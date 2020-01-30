%% Generates forcing for each water year using TopoSCALE.
close all; clear all; clc;
load('../DEM/Lakes_Terrain_Parameters.mat');
load('all_NLDAS.mat'); 


dataset='NLDAS';

yrs=2016:2019;
ny=numel(yrs);

for yr=1:ny
    yris=yrs(yr);
    fprintf('\n Performing the TopoSCALE algorithm on %s data for WY %d \n',dataset,yris);
    t0=datenum(sprintf('30-Sep-%d',yris-1));
    tN=datenum(sprintf('02-Oct-%d',yris));
    dt=1/24;
    t=t0:dt:tN; 
    f=TopoSCALE_NLDAS(t,tp,nldas);
    target=sprintf('TS_%s_%d.mat',dataset,yris);
    save(target,'f','-v7.3');
end
    
