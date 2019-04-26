function [xk1,uk] = trueModel(xk, uk, Ts)
% -------------------------------------------------------------------------
% Simulating the continuous time system
% -------------------------------------------------------------------------
% Inputs:
%   xk: states at time k
%   uk: input at time k
%   Ts: sampling time
%
% Outputs:
%   xk1: states at time k+1
%   uk: input at time k
%
[t,xk] = simulate_cont(uk,xk,[0 Ts]);
xk1 = xk(end,:)';
% choose other ODE solvers such as ode23, ode45 for better accuracy or
% ode15s and ode23s for stiff ODEs
end

function [f] = mgcm_model(x,u)
% -------------------------------------------------------------------------
% Discrete-time nonlinear Moore-Greitzer Compressor dynamic model
% -------------------------------------------------------------------------
% 4 states (x): 
%   x1 - mass flow
%   x2 - pressure rise
%   x3 - throttle
%   x4 - throttle rate
%
% 1 inputs: (u)
%
wn=sqrt(1000); % resonant frequency
zeta=1/sqrt(2); % damping coefficient
beta=1; % constant >0
x2_c=0; % pressure constant
% Continous time state-space model of the Moore-Greitzer compressor model
f = x;
f(1) = -x(2)+x2_c+1+3*(x(1)/2)-(x(1)^3/2); % mass flow rate
f(2) = (x(1)+1-x(3)*sqrt(x(2)))/(beta^2); % pressure rise rate
f(3) = x(4); % throttle opening rate
f(4) = -wn^2*x(3)-2*zeta*wn*x(4)+wn^2*u; % throttle opening acceleration
end

% -------------------------------------------------------------------------
% Simulating the continuous time system /w ODE solvers
% -------------------------------------------------------------------------
function [t, y] = simulate_cont(u, x0, tt)
    [t, y] = ode23(@(t, x)rhs(t, x, u), tt, x0);
end

function f = rhs(t, x, u)
    f = mgcm_model(x,u);
end

