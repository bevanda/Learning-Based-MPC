function [xk1] = trueModel(xk, uk)
%% Discrete-time linear dynamic model 
% 
% 4 states (x): 
%
% 1 inputs: (u)
% -------------------------------------------------------------------------
% Simulating the continuous time system
% -------------------------------------------------------------------------
% Repeat application of Euler method sampled at Ts/M.
Ts = 0.01;
M = 10;
delta = Ts/M;
xk1 = xk;
for ct=1:M
    xk1 = xk1 + delta*dynamics(xk1,uk);
end
% Note that we choose the Euler method (first oder Runge-Kutta method)
% because it is more efficient for plant with non-stiff ODEs.  You can
% choose other ODE solvers such as ode23, ode45 for better accuracy or
% ode15s and ode23s for stiff ODEs.  Those solvers are available from
% MATLAB.
end

function [f] = dynamics(x, u)
%% Continuous-time nonlinear dynamic model of a pendulum on a cart
%
% 4 states (x): 
% 
% 
% 
% 
% 
% 1 inputs: (u)
%   
%
%
% dxdt is the derivative of the states.
% [A B C D] are state space matrices linearized at the current operating point.
%
wn=sqrt(1000); % resonant frequency
zeta=1/sqrt(2); % damping coefficient
beta=1; % constant >0
x2_c=0; % pressure constant

%% Continous time state-space model of the Moore-Greitzer compressor model
f = x;
f(1) = -x(2)+x2_c+1+3*(x(1)/2)-(x(1)^3/2); % mass flow rate
f(2) = (x(1)+1-x(3)*sqrt(x(2)))/(beta^2);  % pressure rise rate
f(3) = x(4);                               % throttle opening rate
f(4) = -wn^2*x(3)-2*zeta*wn*x(4)+wn^2*u;   % throttle opening acceleration

end