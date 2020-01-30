%% main.m
% Main script for running the ensemble-based snow reanalysis for 4 water
% years over the Mammoth Lakes basin.
close all; clear all;
load('../DEM/Mammoth_Lakes_NLCD.mat');
load('../DEM/Lakes_Terrain_Parameters.mat');
slp=tp.slp(tp.mask);
svf=tp.svf(tp.mask);

flim=0.25;
CF=NLCD.grid1.CF./1e2; % Multiply by zero to turn off canopy normalization.
WF=NLCD.grid1.WF;
assimall=1; % Flag to assimilate all retrievals (1=yes, 0=no).
yrs=2016:2019; ny=numel(yrs);
R=0.1.^2;
p.Ne=1e2;

for yr=1:ny
    yris=yrs(yr);
    load(sprintf('../Forcing/tagg_forcing_%d.mat',yris));
    f=fa;
    if assimall==1
        load(sprintf('../Retrievals/Gridded_Obs_All/obs_stack_%d.mat',yris));
    else
        load(sprintf('../Retrievals/Gridded_Obs_noS2/obs_stack_%d.mat',yris));
    end
    CF=CF(f.mask);
    CF=repmat(CF,1,size(obs.fSCA,2));
    obs.fSCA=obs.fSCA./(1-CF); obs.fSCA(obs.fSCA>1)=1; % Canopy normalization for fSCA retrievals.
    p.yearis=yris;
    setpars;
    Np=sum(f.mask(:));
    WFp=WF(f.mask);
    
    % Define the output structure.
    tmp=nan(Np,numel(p.t),3,'single');
    r.pbs.fSCA=tmp; r.pbs.D=tmp;
    r.prior.fSCA=tmp; r.prior.D=tmp;
    r.es.fSCA=tmp; r.es.D=tmp;
    r.esmda.fSCA=tmp; r.esmda.D=tmp;
    
    %% Loop over grid cells.
    decade=(0.1:0.1:1).*1e1;
    for i=1:1%Np
        disp(i./Np);%
        p.slp=slp(i);
        p.svf=svf(i);
        rsc=round(1e1.*i/Np);
        here=(decade-rsc)==0;
        if any(here)
            fprintf('\n %d fraction complete=%4.2f \n',yris,rsc./1e1);
            decade(here)=NaN;
        end
        % Prepare forcing for this grid cell.
        f.t=fa.t;
        these=f.t>=p.t(1)&f.t<=p.t(end);
        f.Qh=fa.Qh(i,these); f.Qe=fa.Qe(i,these);
        f.SW=fa.SW(i,these); f.LW=fa.LW(i,these);
        f.Ps=fa.Ps(i,these); f.Pr=fa.Pr(i,these);
        f.Ta=fa.Ta(i,these);
        o.R=R;
        o.fSCA=obs.fSCA(i,:); o.t=obs.t(:);
        if ~assimall
            o.fSCA(obs.type~=2)=NaN;
        end
        these=~isnan(o.fSCA);
        o.fSCA=o.fSCA(these); o.t=o.t(these);
        
        try
            if WFp(i)>flim
                disp('Wrong land cover');
            else
                % Run the reanalysis for this grid cell.
                [ana,p]=EnSnoRe(o,p,f);
                % PBS
                r.pbs.fSCA(i,:,:)=ana(1).fSCA;
                r.pbs.D(i,:,:)=ana(1).D;
                % Prior
                r.prior.fSCA(i,:,:)=ana(2).fSCA;
                r.prior.D(i,:,:)=ana(2).D;
                % ES
                r.es.fSCA(i,:,:)=ana(3).fSCA;
                r.es.D(i,:,:)=ana(3).D;
                % ES-MDA
                r.esmda.fSCA(i,:,:)=ana(4).fSCA;
                r.esmda.D(i,:,:)=ana(4).D;
                
            end
        catch
            error('Crashed at i=%d',i);
        end
    end
    
    r.x=f.x; r.y=f.y; r.mask=f.mask; r.t=p.t;
    save(sprintf('result_%d.mat',yris),'r','-v7.3');

end



    
 
