function [Ms, Ps, Pr] = getPMrates(f,c,alb,n,p,bs) 
    %% Net accumulation and ablation at the current time step.
    Ps=f.Ps(n)./1e3; Pr=f.Pr(n)./1e3; % Convert to m/day
    Ta=f.Ta(n);
    if any(Ta<0) % Ensure that temperature is in Kelvin.
        Ta=Ta+c.Tm;
    end
    Qe=f.Qe(n); Qh=f.Qh(n);
    SWDnet=(1-alb).*f.SW(n)+(1-alb).*alb.*0.5.*(1+cos(p.slp)).*(1-p.svf).*f.SW(n); % Second term is estimate of radiation reflected by terrain.
    LWDnet=f.LW(n)-c.sbc.*c.eps_snow.*(c.Tm.^4); % Net longwave.
    Qp=(1/c.spd).*(bs.*Ps.*c.c_i.*c.rho_w.*min(Ta-c.Tm,0)+...
        bs.*Pr.*(c.Lf.*c.rho_w+c.c_w.*c.rho_w.*max(Ta-c.Tm,0))); % Heat advected by precipitation.
    QT=(Qh+Qe); % Net turbulent flux. 
    QR=LWDnet+SWDnet; % Net radiation.
    RES=QR-QT+Qp; % Energy balance residual.
    Ms=RES.*(c.spd.*p.dt./(c.Lf.*c.rho_w)); % Melt rate. 
    Es=Qe.*c.spd.*p.dt/(c.Ls.*c.rho_w); 
    Ms=Ms+Es; 
    
end
