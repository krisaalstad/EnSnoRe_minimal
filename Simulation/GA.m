function [ Xt ] = GA( X , bounds , delta, F )
%% GA: Quasi bijective analytical Guassian anamorphosis functions used in the 
%ESSDA framework to improve the optimality of the KF analysis equation by
%transforming relevant variables to an unbounded space. The F status 
%variable (0/1) specifies whether or not this is a forwards transform to
%Gaussian space or inverse transform to physical space.

%%
%{ 
%-------------------------Description-------------------------------------
Abreviations in the description: 
              Ns = Number of state variables  (rows).
              Ne = Number of ensemble members (columns).
              L = Finite non-negative lower bound  (0<=L<Inf).
              U = Finite non-negative uppper bound (0<=L<U<Inf).
              

Inputs: 

       X = Ns x Ne ensemble state matrix with entries to be transformed.

       bounds = Ns x 2 matrix specifying the bounds of the state in
                physical space. Respective variables/parameter bounds are 
                given in each row. Cases:
            
                    1) Unbounded   [-Inf Inf]: No transformation.
                        Example variable: Velocity components.
            
                    2) Lower bounded  [L Inf]: Scaled log transform.
                        Example variable: Tracer concentration.
            
                    3) Double bounded [L   U]: Scaled logit transform.
                        Example variable: Snow cover fraction.

       delta = Fudge factor (0<delta<<1) to avoid discontinuities at the 
                boundaries in the scaled logit transform. 
                Recommend value is 1e-2.

        F = Transform status. F=1 performs a transform from the
        physical space, F=0 performs a transform to the physical space.

Output:
 
        Xt = Ns x Ne transformed ensemble state matrix.

%--------------------------------------------------------------------------
%}
%%
%{
%-----------------------Example--------------------------------------------
%% Inflation in transformed space.
clear all; close all;
Ne=1e6;
x=repmat([0.0;0.5;1.00],1,Ne);  % Dummy bounded ensemble state matrix. 
Ns=size(x,1);                   % Number of state variables.  
delta=1e-2;                     % Fudge factor.
bounds=repmat([0 1],Ns,1);      % Bounds of the dummy matrix (same for both
                                % variables in this case).


% Forwards transform to unbounded space and inflation.
xts=GA(x,bounds,delta,1);       % Unbounded ensemble state matrix.
sdp=0.5;                        % Standard deviation in unbounded space.
xtsp=xts+sdp*randn(size(xts));

% Inverse transform of the inflated and unbounded state to
% bounded physical space.
xp=GA(xtsp,bounds,delta,0);

% Visualization of pdfs normalized by their maxima.
f1=figure(1);
set(f1,'Units','Normalized','OuterPosition',[0 0 1 1]);
cs=0.5.*eye(Ns);
bins=0:0.01:1;
for n=1:Ns
    ht=hist(xp(n,:),bins,'hist');
    b(n)=bar(bins,ht./max(ht));
    set(b(n),'BarWidth',1,'FaceColor',cs(n,:),'EdgeColor',cs(n,:));
    alpha(0.2);
    hold on;
end
leg=legend(b(:),{'$x_1$','$x_2$','$x_3$'},'Location','NorthOutside',...
    'Orientation','Horizontal','Interpreter','Latex','FontSize',24);
xlabel('$x$','Interpreter','Latex','FontSize',24);
ylabel('$p_x/max(p_x)$','Interpreter','Latex','FontSize',24);
xlim([min(min(xp)) max(max(xp))]);
box on; grid on;
set(gca,'TickDir','out');

%--------------------------------------------------------------------------
%}


%%
%-----------------------Main routine---------------------------------------

Xt=zeros(size(X)); % Transformed ensemble state matrix pre-allocation. 
Ns=size(X,1);      % Number of state variables.

for n=1:Ns % Loop over state variables, the ensemble loop is vectorized.
    if all(isinf(bounds(n,:)))     % Case 1) Unbounded state variable.  
        Xt(n,:)=X(n,:);
    elseif any(isinf(bounds(n,:))) % Case 2) Lower bounded state variable.
        Xt(n,:)=slog(X(n,:),bounds(n,1));
    else                           % Case 3) Double bounded state variable. 
        Xt(n,:)=slogit(X(n,:),bounds(n,1),bounds(n,2));
    end    
end


%-----------------------Nested functions-----------------------------------

%% Log transforms.
% Forwards transform (F=1): maps x in [L,Inf] to xts [-Inf,Inf].
% Inverse transform (F=0): maps x in [-Inf,Inf] to xts [L,Inf].
    function xts=slog(x,L)
	l=L-delta; 
        if F
            x(x<L)=L;
            xs=x-l;
            xts=log(xs);
        else
            xt=exp(x);
            xts=xt+l;
            xts(xts<L)=L;
        end
    end


%% Logit transforms.
% Forwards transform (F=1): maps x in [L,U] to xts [-Inf,Inf].
% Inverse transform (F=0): maps x in [-Inf,Inf] to xts in [L,U].
    function xts=slogit(x,L,U)
        l=L-delta; u=U+delta;
        if F
            xs=(x-l)./(u-l);
            xts=log(xs)-log(1-xs);
        else
            xt=(u-l)./(1+exp(-x))+l;
            xts=xt; xts(xts<L)=L; xts(xts>U)=U;
        end
    end
end
