function [Qe,Qh]=MO_turbulent_fluxes(c,T,q,U,ps)
%% Monin-Obukhov turbulent heat flux calculations.
% Simplified scheme that assumes the surface is at 0^degC.
% Based on Westermann et al. 2016.
% See also Paulson 1970 for more details https://doi.org/10.1175/1520-0450(1970)009%3C0857:TMROWS%3E2.0.CO;2


[Qe,Qh]=deal(zeros(size(T))); % Pre allocation.
Nt=size(Qh,2); % Number of time steps.
Lstar0=1e3; % Initial Obukhov length (not that important).
rho=ps./(T.*c.Rd); % Air density.
z0=c.z0_snow;
% Time loop.
for n=1:Nt
    disp(n./Nt);
    if n==1
        Lstar=Lstar0;
    else
        ustar=U(:,n).*c.Kappa./...
            (log(c.zU/z0)-psi_m(c.zU./Lstar,z0./Lstar));
        Lstar=rho(:,n).*c.cp.*T(:,n).*(ustar.^3)./...
            (c.Kappa.*c.g.*(Qh(:,n-1)+c.evd.*c.cp.*T(:,n).*Qe(:,n-1)./c.Ls));
    end
    ustar=U(:,n).*c.Kappa./...
            (log(c.zU/z0)-psi_m(c.zU./Lstar,z0./Lstar));
    ra=((c.Kappa.*ustar).^(-1)).*(log(c.zTq./z0)-psi_h(c.zTq./Lstar,z0./Lstar));
    if any(ra<0)
        error('error aerodynamic resistance should be positive');
    end
    Qh(:,n)=rho(:,n).*c.cp.*(c.Tm-T(:,n))./ra;
    qsis=qsat0(c,ps(:,n));
    Qe(:,n)=rho(:,n).*c.Ls.*(qsis-q(:,n))./ra;
end


    function qs=qsat0(c,ps)
        %% Saturation Specific humidity at 0 degrees Celsius
        % as a function of pressure [in Pa].
        es=c.e0; % Saturation vapor pressure at 0^\degC. 
        mrvs=c.evd.*es./(ps-es); % Mixing ratio at saturation.
        qs=mrvs./(1+mrvs); % Saturation specific humidity.
        % N.B. qs is the same w.r.t. both ice and water at 0^\degC.
    end



    function res=psi_h(zeta1,zeta2)
        %% Integrated stability function for heat.
        % Follows notation in Westermann et al. 2016 integrating FROM zeta2
        % to zeta1, where zeta1=z/L_* and zeta2=z0/L_*
        stable=zeta1>=0; % Stable?
        
        % Note that psi_H is the integral of (1-phi_H)/zeta from zeta2 to
        % zeta1 where phi_H is the nondimensional gradient of potential
        % temperature in the surface layer.
        
        %% Unstable case (zeta<0)
        % Using psi_H=(1-12*zeta)^(-1/2) from page 77 in Högström (1988) in
        % the unstable case. Use of this nondimensional gradient avoids a
        % discontinuity between the stable and unstable regime at zeta=0
        % (i.e. neutral stability).
        psi_unstable=-2.*(log(abs((1-12.*zeta1).^(-1/2)))-log(abs((1-12.*zeta1).^(-1/2)+1)))...
            +2.*(log(abs((1-12.*zeta2).^(-1/2)))-log(abs((1-12.*zeta2).^(-1/2)+1)));
        
        %% Stable case (zeta>0)
        % Following Grachev et al. (2007) in the stable case.
        ah=5;
        bh=5;
        ch=3;
        Bh=sqrt(ch^2-4);
        psi_stable=(-bh/2).*log(1+ch.*zeta1+zeta1.^2)+((-ah/Bh)+bh*ch/(2*Bh)).*(log((2.*zeta1+ch-Bh)./(2.*zeta1+ch+Bh))-log((ch-Bh)/(ch+Bh)))...
            -(-bh/2).*log(1+ch.*zeta2+zeta2.^2)-((-ah/Bh)+bh*ch/(2*Bh)).*(log((2.*zeta2+ch-Bh)./(2.*zeta2+ch+Bh))-log((ch-Bh)/(ch+Bh)));
        
        
        %% Output
        res=(~stable).*psi_unstable+(stable).*psi_stable; % result
    end


    function res=psi_m(zeta1,zeta2)
        %% Integrated stability function for momentum.
        % Follows notation in Westermann et al. (2016) integrating FROM zeta2
        % to zeta1, where zeta1=z/L_* and zeta2=z0/L_*
        stable=zeta1>0; % Stable?
        
        % Note that psi_m is the integral of (1-phi_m)/zeta from zeta2 to
        % zeta1 where phi_m is the nondimensional gradient of wind speed
        % in the surface layer.
        
        
        %% Unstable case.
        % Using psi_m=(1-19*zeta)^(-1/4) from page 73 in Högström (1988) in
        % the unstable case.
        psi_unstable=-4.*(0.5.*(-0.5.*log(abs((1-19.*zeta1).^(-1/2)+1))-atan((1-19.*zeta1).^(-1/4)))+log(abs((1-19.*zeta1).^(-1/4)))-0.5.*log(abs((1-19.*zeta1).^(-1/4)+1)))...
            +4.*(0.5.*(-0.5.*log(abs((1-19.*zeta2).^(-1/2)+1))-atan((1-19.*zeta2).^(-1/4)))+log(abs((1-19.*zeta2).^(-1/4)))-0.5.*log(abs((1-19.*zeta2).^(-1/4)+1)));
        
        
        %% Stable case.
        % Following Grachev et al. (2007) in the stable case.
        am=5;
        bm=am/6.5;
        Bm=((1-bm)/bm)^(1/3);
        x1=(1+zeta1).^(1/3);
        x2=(1+zeta2).^(1/3);
        psi_stable=(-3*am/bm).*(x1-1)+(am*Bm/(2*bm)).*(2*log((x1+Bm)./(1+Bm))-log((x1.^2-x1.*Bm+Bm^2)./(1-Bm+Bm^2))+sqrt(12).*(atan((2.*x1-Bm)./(sqrt(3)*Bm))-atan((2-Bm)/(sqrt(3)*Bm))))...
            -(-3*am/bm).*(x2-1)-(am*Bm/(2*bm)).*(2*log((x2+Bm)./(1+Bm))-log((x2.^2-x2.*Bm+Bm^2)./(1-Bm+Bm^2))+sqrt(12).*(atan((2.*x2-Bm)./(sqrt(3)*Bm))-atan((2-Bm)/(sqrt(3)*Bm))));
        
        
        %% Output.
        res=(~stable).*psi_unstable+stable.*psi_stable;
        
    end








end