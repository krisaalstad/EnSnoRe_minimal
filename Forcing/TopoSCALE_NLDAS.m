function fout=TopoSCALE_NLDAS(t,tp,fin)
%% Apply the TopoSCALE routine to downscale NLDAS data.
% Based on the TopoSCALE algorithm described in Fiddes and Gruber (2014). 
% Code by K. Aalstad (June 2019).

% Scale factors:
fout.wind_sf=1e-2; 
fout.q_sf=1e-6;
fout.ps_sf=1e2;
fout.rad_sf=1e-1;
fout.T_sf=1e-2;
fout.P_sf=1e-2;

Nt=numel(t);
fout.t=t;
fout.x=tp.x;
fout.y=tp.y;
fout.mask=tp.mask;
fout.utmz=tp.utmz;

%% Canopy parameters. 
% Following Bair 2016 for coniferous trees.
hc=16; % Height of canopy (m).
tauc=0.3; % Canopy transmissivity.
muc=0.033; % Canopy extinction coefficient. 

%% To UTM
[LAT,LON]=meshgrid(fin.lat,fin.lon);
utmstruct=defaultm('utm');
utmstruct.zone=tp.utmz;
utmstruct.geoid=wgs84Ellipsoid;
utmstruct=defaultm(utmstruct);
[Xc,Yc]=mfwdtran(utmstruct,LAT,LON);

[X,Y]=meshgrid(tp.x,tp.y);
X=X(tp.mask); Y=Y(tp.mask);

%% Output structure
temp=zeros(numel(X),Nt,'single');
[fout.T,fout.u,fout.v,fout.q,fout.P,fout.LW,fout.SW,fout.ps]=deal(temp);
Zc=temp;


%% Compute the IDW weights
Np=numel(X); w=zeros(Np,numel(Xc));
for k=1:Np
    xis=X(k); yis=Y(k);
    % Distances
    d=sqrt((xis-Xc).^2+(yis-Yc).^2);
    dlim=sort(d(:)); dlim=dlim(4);
    % Only include 4 nearest.
    wis=1./(d.^2).*(d<=dlim);
    wis=wis./sum(wis(:));
    w(k,:)=wis(:);
end

%% Terrain parameters
slp=tp.slp(tp.mask);
asp=tp.asp(tp.mask);
svf=tp.svf(tp.mask);
nb=numel(tp.hbins); h=zeros(numel(X),nb);
for b=1:nb
    htemp=tp.h(:,:,b); 
    h(:,b)=htemp(tp.mask);
end


%% Time loop.
outstatus=1:1:100;
for n=1:Nt
    
    rpc=round(1e2*(n/Nt));
    here=outstatus==rpc;
    if any(here)
        outstatus(here)=NaN;
        fprintf('\n %d percent of downscaling complete \n',rpc);
    end
    
    here=fin.t==t(n);
    here=find(here==1); % Faster to use index than boolean array.
    if ~any(here)
        error('Missing NLDAS data for this time period');
    end

    %% Disaggregate to high res grid using IDW
    % Vectorized :)
    tmp=fin.z;
    Zc=w*tmp(:); % Disaggregated NLDAS elevations.
    
    tmp=fin.P(:,:,here);
    fout.P(:,n)=w*tmp(:); 
    
    tmp=fin.T(:,:,here);
    fout.T(:,n)=w*tmp(:);
    
    tmp=fin.u(:,:,here);
    fout.u(:,n)=w*tmp(:);
    
    tmp=fin.v(:,:,here);
    fout.v(:,n)=w*tmp(:);
    
    tmp=fin.q(:,:,here);
    fout.q(:,n)=w*tmp(:);
    
    tmp=fin.SW(:,:,here);
    fout.SW(:,n)=w*tmp(:);
    
    tmp=fin.LW(:,:,here);
    fout.LW(:,n)=w*tmp(:);
    
    tmp=fin.p(:,:,here);
    fout.ps(:,n)=w*tmp(:);
    
    %% Temperature.
    % Use fixed-lapse rate following Girotto et al. (2014)
    dz=tp.z(tp.mask)-Zc;
    Tc=fout.T(:,n); % Save old value (prior to lapse rate correction).
    fout.T(:,n)=fout.T(:,n)-5.*dz./1e3; % - 5 deg / km.
    
    %% Pressure.
    Tbar=(Tc+fout.T(:,n))/2;
    R=287; g=9.81;
    H=R.*Tbar./g;
    pc=fout.ps(:,n);
    fout.ps(:,n)=pc.*exp(-dz./H); % Hypsometric equation.
    
    %% Humidity.
    % Use lapse rate following Liston and Elder (2006)
    vp=fout.q(:,n).*fout.ps(:,n)./0.622;
    vpc=vp;
    a1=7.625;b1=243.04;c1=610.94;
    Td=b1.*log(vp./c1)./(a1-log(vp./c1)); % From Lawrence (2005) in BAMS.
    Td=Td-(0.4.*240/17.502).*dz./1e3; % Apply lapse rate correction, numbers from Liston and Elder.
    Tf=fout.T(:,n)-273.15;
    c2=a1.*Tf./(b1+Tf);
    RH=exp((a1.*Td-c2.*Td-b1.*c2)./(b1+Td)); % As a fraction.
    vps=c1.*exp(a1.*Tf./(b1+Tf));
    vp=RH.*vps;
    wmr=0.622.*vp./(fout.ps(:,n)-vp);
    fout.q(:,n)=wmr./(1+wmr);
    
    %% Longwave
    
    % Use the vapor pressure and temperature to calculate clear sky
    % emssivity at grid and subgrid. 
    x1=0.43; x2=5.7;
    cef=0.23+x1.*(vp./fout.T(:,n)).^(1/x2); 
    cec=0.23+x1.*(vpc./Tc).^(1/x2);
    
    % Diagnose the all sky emissivity at grid.
    sbc=5.67e-8; 
    aec=fout.LW(:,n)./(sbc.*Tc.^4);
    
    % Calculate the "cloud" emissivity at grid, assume this is the same at
    % subgrid.
    deltae=aec-cec;
    
    % Use the former cloud emissivity to compute the all sky emissivity at subgrid.
    aef=cef+deltae;
    LWf=aef.*sbc.*fout.T(:,n).^4;
    
    % Scale LW with terrain configuration, considering both occlusion by
    % and emissions from the surrounding terrain. From Dozier & Frew 1990.
    LWf=svf.*LWf+0.5.*(1+cos(slp)).*(1-svf).*0.99.*sbc.*(273.15.^4); 
   
    % Apply canopy correction.
    LWf=tp.cc.*(LWf.*tauc+sbc.*0.96.*fout.T(:,n).^4)+(1-tp.cc).*LWf;
    
    % Set longwave.
    fout.LW(:,n)=LWf;
    
    %% Shortwave.
    SWc=fout.SW(:,n);
    
    there=tp.t==t(n);
    sz=tp.szen(there);
    sz=sz.*(sz<pi/2); % Sun might be below the horizon. Solar zenith of 90 is below the horizon.
    muz=cos(sz); 
    S0=1362; % W/m^2
    SWtoa=S0*cos(tp.szen(there)); % Incoming SW at TOA. Ignoring eccentricity in orbit for simplicity.
    SWtoa(SWtoa<0)=0;
    kt=SWc./SWtoa;
    kt(isnan(kt))=1; % division by zero (night) will lead to NaN.
    
    % Calculate the diffuse fraction following the regression of Ruiz-Arias
    kd=0.952-1.041.*exp(-1.*exp(2.3-4.702.*kt)); 
    kd=kd.*(kd>0);
    
    % Use this to calculate the downwelling diffuse and direct shortwave radiation at grid.
    SWcdiff=kd.*SWc;
    SWcdir=(1-kd).*SWc;
    
    % Use the above with the sky-view fraction to calculate the downwelling diffuse
    % shortwave radiation at subgrid. 
    SWfdiff=svf.*SWcdiff;
    
    % Calculate the broadband absorption coefficient.
    ka=(g.*muz./(pc)).*log(SWtoa./SWcdir);
    ka(isnan(ka)|ka<0)=0; 
    
    % Find the direct component at subgrid. 
    SWfdir=SWtoa.*exp(-ka.*fout.ps(:,n)./(g*muz)); 
    
    % Illumination angles.
    saz=tp.saz(there); % solar azimuth.
    cosis=muz.*cos(slp)+sin(sz).*sin(slp).*cos(saz-asp); % cosine of illumination angle at subgrid.
    cosis(cosis<0)=0; % Self-shadowing |saz-tp.asp|>90
    cosic=muz; % cosine of illumination angle at grid (assumes slope=0).
    
    % Shadow mask
    these=abs(tp.hbins-saz)==min(abs(tp.hbins-saz));
    temp=1:nb;
    nnbin=min(temp(these));
    sel=max(rad2deg(pi/2-sz),0);
    shade=h(:,nnbin)>sel;
    SWfdir=SWfdir.*(cosis./cosic).*(1-shade); % Terrain corrected direct shortwave at subgrid.
    
    % Canopy corrections.
    SWfdir=tp.cc.*(SWfdir.*exp(-muc.*hc/cos(sz)))+(1-tp.cc).*SWfdir;
    SWfdiff=tp.cc.*(SWfdiff.*tauc)+(1-tp.cc).*SWfdiff;
    
    % Compute total incoming shortwave at subgrid.
    fout.SW(:,n)=SWfdir+SWfdiff;
    
    
end


%% Compress the results.
fout.u=int16(fout.u./fout.wind_sf);
fout.v=int16(fout.v./fout.wind_sf);
fout.q=uint16(fout.q./fout.q_sf);
fout.ps=uint16(fout.ps./fout.ps_sf);
fout.SW=uint16(fout.SW./fout.rad_sf);
fout.LW=uint16(fout.LW./fout.rad_sf);
fout.T=int16((fout.T-273.15)./fout.T_sf);
fout.P=uint16(fout.P./fout.P_sf);


end







