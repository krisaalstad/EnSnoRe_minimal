function [mu_new,D_m_new,fSCA,D_a,albs_new]=...
    SSM(mu_cur,D_m_cur,albs_cur,NA,chi,p,albsmin)
%% SSM: Simple Snow Model
%   Runs one timestep of SSM. The model is fully
%   vectorized allowing it to be run for an ensemble.
%   Inputs (_cur = value at current timestep):
%       mu_cur   = Accumulated peak mean snow depth (m).
%       D_m_cur  = Accumulated melt depth (m).
%       albs_cur= Snow albedo.
%       NA      = Net accumulation (m/day). 
%       chi     = Coefficient of variation of the snow distribution.
%       p       = Structure of parameters (defined in setpars.m)
%       albsmin = Minimum snow albedo. 
%   Outputs (_new = value at next time step).
%       mu_new    = Accumulated peak mean snow depth (m).
%       D_m_new   = Accumulated melt depth (m).
%       fSCA      = Fractional snow-covered area (-).
%       D_a       = Mean SWE depth (m).
%       albs_new  = Snow albedo.

    
    %% Update accumulated melt depth (m).
    issnow=mu_cur>0;
    D_m_new=max(D_m_cur-NA,0).*issnow; % Important to include, otherwise you will carry a huge melt depth from the summer.
    
    
    %% Update peak accumulated depth (m).
    mu_new=(mu_cur+(max(NA-(D_m_new),0))); % Avoids double counting.
    issnow=mu_new>0;

     %% Update fSCA (-) and average SWE depth (m).
     if any(issnow)
        [fSCA,D_a]=fSCA_scheme();
     else
         [fSCA,D_a]=deal(zeros(size(NA))); 
     end
     fSCA(issnow==0)=0; D_a(issnow==0)=0;
    
    %% Update the snow albedo.
    albs_new=SALB_scheme();
    
    % Resetting is vital, mu and Dm are just accounting variables.
    issnow=fSCA>0;
    D_m_new=D_m_new.*issnow;
    mu_new=mu_new.*issnow;
    
    %% fSCA scheme.
    function [fSCA,D_a]=fSCA_scheme()
        % Based on Liston (2004; JClim).
        % Calculates the fSCA given an accumulated melt depth increment 
        % D_m_new-D_m_cur and an initial (premelt) lognormal snow distribution 
        % with mean mu_new and coefficient of variation (CV=standard dev/mean).
        
        sdt=sqrt(log(1+chi.^2));                       % SD of the normal distribution of ln(D).
        mut=log(mu_new)-0.5.*sdt.^2;	               % Mean of the normal distribution of ln(D).
        mut(mut==-Inf)=0;			                   % Needed in the case that D_m=Inf as Inf-Inf=NaN.
        zDmnew=(log(D_m_new)-mut)./(sqrt(2).*sdt);     % Shorthand variable given current melt rate.
        
        fSCA=0.5.*erfc(zDmnew);                        % Diagnose fSCA.
        
        fSCA(fSCA<1e-2)=0;                             % Set to zero if <0.01.
        
        D_a=0.5.*exp(mut+0.5.*sdt.^2).*...             % Diagnose grid cell mean SWE.
            erfc((zDmnew-sdt)./sqrt(2))...
            -fSCA.*D_m_new;
        
        D_a(isnan(D_a))=0; D_a(D_a<0)=0;
        D_a(fSCA==0)=0;
    end

    %% Snow albedo scheme.
    function albs_new=SALB_scheme()
        % Follows equation A7 in Dutra et al. (2010);
        albs_new=(NA>0).*(albs_cur+min(1,(NA.*p.dt./p.taus)).*(p.albs_max-albs_cur))+...
            (NA==0).*(albs_cur-p.taua.*p.dt)+...
            (NA<0).*((albs_cur-albsmin).*exp(-p.tauf.*p.dt)+albsmin);
        albs_new(albs_new<albsmin)=albsmin(albs_new<albsmin);
    end
     
end



