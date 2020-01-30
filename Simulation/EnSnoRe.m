function [ana,p]=EnSnoRe(o,p,f)
%% [ana,p]=EnSnoRe(observations,parameters,forcing)
% EnSnoRe: Ensemble-based Snow Reanalysis routine.
% Runs the reanalysis for a given water year and model grid cell.
% By K. Aalstad (Dec 2016, last revised Jan 2020).
% Requires:
% 1) An observation structure 'o' with the following fields specified:
%                   o.fSCA  = 1D Array of fSCA observations.
%                   o.R     = Observation error covariance matrix.
%                   o.t     = 1D Array of timestamps for these
%                                   observations.
% 2) A parameter structure 'p' with the following fields specified:
%                   pars.yearis = water year.
% Remaining parameters should be specified in setpars.m which is run
% internally in a call to EnSnoRe.
% 3) A forcing structure 'f' with the following fields specified:
%                   f.Qh      = Sensible heat flux [Wm^-2] 
%                   f.Qe      = Latent heat flux [Wm^-2]
%                   f.SW      = Shortwave radiation flux [Wm^-2]. 
%                   f.LW      = Longwave radiation flux [Wm^-2]. 
%                   f.Ps      = Snowfall rate [mm/day]. 
%                   f.Pr      = Rainfall rate [mm/day]. 
%                   f.Ta      = Air temperature [K].
% All the forcing fields are expected to be time series (1D arrays) for the
% grid cell in question.


%% Setup the experiment.

% Set the remaining parameters.
setpars;

% Get constants.
get_constants;

% Initialize the seed.
rdraw=round(rand*1e5); 
rng(rdraw); 

% Allocate the observations.
fSCAo=o.fSCA; to=o.t;                    
p.R=o.R; p.No=numel(fSCAo); clear obs;


%% Define the output array structure
nm=4;
methods={'PBS','Prior','ES','ES-MDA'};
for m=1:nm
    ana(m).method=methods{m};
end

%% Pre-Allocation.
[ana(:).chi,ana(:).bp,ana(:).bm,ana(:).mu,ana(:).albsmin]=deal(zeros(1,p.Ne));
[ana(:).fSCA,ana(:).D]=deal(zeros(p.Nt,numel(p.pct)));
[pfSCA,pD]=deal(zeros(p.Nt,p.Ne));   
HX=zeros(p.No,p.Ne);
y=deal(zeros(p.No,1));
stopit=0;
l=0; % Iteration counter.



%% Assimilation cycles.
while ~stopit
    
    
    % Define the number of independent parameter ensembles going into the
    % given forwards run (for iteration l):
    npe=(l~=1)+(l==1).*2; % 2 for the second iteration (ES, ES-MDA), 1 otherwise (EITHER prior [l==0] or ES-MDA [l>1 & l~=2]).
    
    no=1; % Observation counter.
    
    %% Prediction loop
    for n=0:p.Nt-1 
        
        %% Initialization.
        if ~n
            
            if ~l
                % Prior parameter ensemble.
                chi=GA(GA(p.mchi,p.chib,0,1)+p.sdtchi.*randn(1,p.Ne),p.chib,0,0);
                bm=GA(GA(p.mbm,p.bmb,0,1)+sqrt(log(1+p.cvbm.^2)).*randn(1,p.Ne),p.bmb,0,0);
                bp=GA(GA(p.mbp,p.bpb,0,1)+sqrt(log(1+p.cvbp.^2)).*randn(1,p.Ne),p.bpb,0,0);
                albsmin=GA(GA(p.malbsmin,p.albsminb,0,1)+p.sdtalbsmin.*randn(1,p.Ne),p.albsminb,0,0);
                for m=1:nm
                    ana(m).chi=chi;
                    ana(m).bp=bp;
                    ana(m).bm=bm;
                    ana(m).albsmin=albsmin;
                end
            elseif l==1
                chi=[ana(end-1).chi ana(end).chi];
                bm=[ana(end-1).bm ana(end).bm];
                bp=[ana(end-1).bp ana(end).bp];
                albsmin=[ana(end-1).albsmin ana(end).albsmin];
            else
                chi=ana(end).chi; bm=ana(end).bm; bp=ana(end).bp; albsmin=ana(end).albsmin;
            end
            
            
            % Initialize state variables as zeros.
            [Dm,mu,albs,fSCA]=deal(zeros(1,npe.*p.Ne)); % Dm=Melt depth, mu=mean pre-melt snow depth,albs =snow albedo.
            
            
            %% Forwards run
        else
            
            %% Forcing
            % Allocate the forcing for the given time step.
            % Diagnose forcing based on SEB and SWB fields.
            alb=fSCA.*albs+(1-fSCA).*0.2; % Effective albedo, implicitly accounts
                                        % for terrain reflection of SW and
                                        % litter/LAI.
            
            % Diagnose precpitation and melt rates.                            
            [Msis, Psis, Pris] = getPMrates(f,c,alb,n+1,p,bp); % Change alb to albs if you dont want to use effective albedo.
            
            % Perturb the forcing.
            Pse=Psis.*bp; Mse=Msis.*bm; Pre=Pris.*bp.*(Dm==0); 
            NAe=((Pse.*(Pse>0)-Mse.*(Mse>0)+Pre.*(Pre>0))).*(p.dt); % Net accumulation in this timestep.
            
            
            %% Model step
            % Run one timestep of SSM
            [mu,Dm,fSCA,D,albs]=SSM(mu,Dm,albs,NAe,chi,p,albsmin);
            
            %% Save outputs.
            if ~l                               
                pfSCA(n+1,:)=fSCA; % Prior fSCA
                pD(n+1,:)=D;       % Prior SWE
                if (n+1)==p.Nt
                    for m=1:2
                        ana(m).mu=mu;
                    end
                end
            end
            if l==1
                % ES percentiles
                for m=3 
                    those=1:p.Ne;
                    ana(m).fSCA(n+1,:)=prctile(fSCA(those),p.pct);
                    ana(m).D(n+1,:)=prctile(D(those),p.pct);
                    if (n+1)==p.Nt
                        ana(m).mu=mu(those);
                    end
                end
                % PBS percentiles   
                sortedD=sortrows([pD(n+1,:)' w'],1);
                cwD=cumsum(sortedD(:,2));
                sortedfSCA=sortrows([pfSCA(n+1,:)' w'],1);
                cwfSCA=cumsum(sortedfSCA(:,2));
                for j=1:numel(p.pct)
                    here=find(abs(cwD-1e-2.*p.pct(j))==min(abs(cwD-1e-2.*p.pct(j))),1);
                    ana(1).D(n+1,j)=sortedD(here,1);
                    here=find(abs(cwfSCA-1e-2.*p.pct(j))==min(abs(cwfSCA-1e-2.*p.pct(j))),1);
                    ana(1).fSCA(n+1,j)=sortedfSCA(here,1);
                end
                % Prior percentiles
                if (n+1)==p.Nt
                    ana(2).D=prctile(pD,p.pct,2);
                    ana(2).fSCA=prctile(pfSCA,p.pct,2);
                    clear pfSCA pD;
                end
            end
            if l==p.Na
                % ES-MDA percentiles
                those=1:p.Ne;
                ana(end).fSCA(n+1,:)=prctile(fSCA(those),p.pct);
                ana(end).D(n+1,:)=prctile(D(those),p.pct);
                if (n+1)==p.Nt
                    ana(end).mu=mu(those);
                end
            end
            
            
            %% Store arrays for the analysis
            % Store observations and predicted observations.
            if no<=p.No
                if p.t(n+1)==to(no)
                    y(no)=fSCAo(no); 
                    HX(no,:)=fSCA(1:p.Ne);
                    no=no+1;
                end
            end
            
        end
    end % End of time loop.
                
    %% Analysis steps.
    if l~=p.Na
        %% ES - MDA
        thetat=GA([ana(end).chi;ana(end).bp;ana(end).bm; ana(end).albsmin],...
            [p.chib; p.bpb; p.bmb; p.albsminb],0,1);
        thetat=fastpEnKF(thetat,HX,y,p.R,p.Na,p.prt_stat);
        theta=GA(thetat,[p.chib; p.bpb; p.bmb; p.albsminb],0,0);
        ana(end).chi=theta(1,:);
        ana(end).bp=theta(2,:);
        ana(end).bm=theta(3,:);
        ana(end).albsmin=theta(4,:);
       
        if ~l
            %% ES
            thetat=GA([ana(end).chi;ana(end).bp;ana(end).bm;ana(end).albsmin],...
                [p.chib; p.bpb; p.bmb; p.albsminb],0,1);
            thetat=fastpEnKF(thetat,HX,y,p.R,1,p.prt_stat);
            theta=GA(thetat,[p.chib; p.bpb; p.bmb; p.albsminb],0,0);
            ana(3).chi=theta(1,:);
            ana(3).bp=theta(2,:); 
            ana(3).bm=theta(3,:); 
            ana(3).albsmin=theta(4,:);
            
            
            %% PBS
            w=PBS(HX,y,p.R); 
            ana(1).w=w;
        end
    else
        stopit=1;
    end
    l=l+1;
end
        
                











