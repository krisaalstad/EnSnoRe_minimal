%% Performs a temporal aggregation of the TopoSCALE forcing. 
% This includes a calculation of the turbulent heat fluxes.

close all; clear all; clc;

yrs=2016:2019; ny=numel(yrs);
dataset='NLDAS';
utcoff=-8;
tfb=1e3;
Train=276.15; Tsnow=272.15;
get_constants;

for yr=1:ny
    yris=yrs(yr);
    fprintf('\n Aggregating forcing WY %d \n',yris);
    tarfile=sprintf('TS_%s_%d',dataset,yris);
    load(tarfile);
    
    
    t=f.t+utcoff./24;
    d=datenum(sprintf('01-Oct-%d',yris-1)):datenum(sprintf('02-Oct-%d',yris));
    Nd=numel(d);
    Np=size(f.T,1);
    
    mask=f.mask; x=f.x; y=f.y; utmz=f.utmz;
    u=double(f.u).*f.wind_sf; v=double(f.v).*f.wind_sf; Ua=sqrt(u.^2+v.^2);
    Ta=double(f.T).*f.T_sf+273.15;
    qa=double(f.q).*f.q_sf;
    ps=double(f.ps).*f.ps_sf;
    Ua(Ua<1)=1;
    [Qe,Qh]=MO_turbulent_fluxes(c,Ta,qa,Ua,ps);
    %Qh(Qh>tfb)=tfb; Qh(Qh<-tfb)=-tfb;
    %Qe(Qe>tfb)=tfb; Qe(Qe<-tfb)=-tfb;
    SW=double(f.SW).*f.rad_sf;
    LW=double(f.LW).*f.rad_sf;
    P=double(f.P).*f.P_sf;
    f_rain=(Ta-Tsnow)./(Train-Tsnow); f_rain(f_rain<0)=0; f_rain(f_rain>1)=1;
    Prain=f_rain.*P; Prain(Prain<0)=0;
    Psnow=P-Prain; Psnow(Psnow<0)=0;
    clear P; clear fout;
    
    fa.x=x; fa.y=y; fa.utmz=utmz; fa.utcoff=utcoff; fa.mask=mask; fa.t=d(1:Nd-1);
    tmp=zeros(Np,Nd-1,'single');
    [fa.Qh,fa.Qe,fa.SW,fa.LW,fa.Ps,fa.Pr,fa.Ta,fa.Ua]=deal(tmp);
    for n=1:Nd-1
        these=(t>=d(n)&t<d(n+1))';
        fa.Qh(:,n)=1.*mean(Qh(:,these),2);
        fa.Qe(:,n)=1.*mean(Qe(:,these),2);
        fa.Ta(:,n)=mean(Ta(:,these),2);
        fa.Ua(:,n)=mean(Ua(:,these),2);
        fa.SW(:,n)=mean(SW(:,these),2);
        fa.LW(:,n)=mean(LW(:,these),2);
        fa.Ps(:,n)=sum(Psnow(:,these),2);
        fa.Pr(:,n)=sum(Prain(:,these),2);
    end
    save(sprintf('tagg_forcing_%d.mat',yris),'fa','-v7.3');
    
end

