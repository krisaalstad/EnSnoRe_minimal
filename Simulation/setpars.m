%% Run parameters for the snow reanalysis framework.

%% Related to time specification. (p.yearis must be specified elsewhere).
p.startdate=datenum(['01-Oct-' sprintf('%d',p.yearis-1)]); % Start of simulation *.
p.enddate=datenum(['01-Oct-' sprintf('%d',p.yearis)]);      % End date of simulation *.
p.dt=1;                                                    % Time step (days).
p.t=datenum(p.startdate):p.dt:datenum(p.enddate);          % Timesteps.
p.Nt=numel(p.t);                                           % Number of time steps.

if numel(fields(p))<30

%% Related to the snow albedo parametrization:
p.albs_max=0.85;                                            % Maximum snow albedo.
p.taus=0.01;                                                % Threshold snowfall for resetting to maximum [m w.e.].
p.taua=0.008;                                               % Time constant for snow albedo change in non-melting conditions [/day].
p.tauf=0.24;                                                % Time constant for snow albedo change in melting conditions [/day].

%% Related to the ensemble:
if ~any(strcmp(fields(p),'Ne'))
    p.Ne=1e2;                                               % Number of ensemble members *.
end
p.pct=[5;50;95];                                            % Ensemble percentiles to select.
p.prt_stat=1;
p.Na=4;


%% Specification of bounds:
p.bpb=[0 Inf];
p.bmb=[0 Inf];
p.chib=[0 1];
p.albsminb=[0.45 0.55]; 



%% Specification of hyperparameters (prior means and variance-related parameters).

% Precipitation bias.
p.cvbp=0.2; % Corresponds to a variance of 0.04. 
p.mbp=1;  % Mean

% Melt bias.
p.cvbm=0.1; % Corresponds to a variance of 0.01.
p.mbm=1; % Mean

% Coefficient of variation.
p.sdtchi=0.4; % Corresponds to a variance of 0.01.
p.mchi=0.4; % Mean

% Minimum snow albedo.
p.sdtalbsmin=1; % Corresponds to a variance of 0.02^2.
p.malbsmin=0.5; % Mean

end
