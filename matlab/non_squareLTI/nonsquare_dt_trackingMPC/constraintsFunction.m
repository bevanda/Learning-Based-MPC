function [c, ceq] = constraintsFunction(u,theta,x,N)
%% Constraint function of nonlinear MPC for pendulum swing-up and balancing control
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%
% Output:
%   c:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints (empty)
%
% Copyright 2016 The MathWorks, Inc.

%% Nonlinear MPC design parameters
c_num = 8;

x_min = -inf; x_max=5.0;
u_min=-inf; u_max=0.3;
%% Inequality constraints calculation
c = zeros(N*c_num,1);
% Apply 2*N cart position constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
uk = u;
for ct=1:N
    % obtain new cart position at next prediction step
    xk1 = getTransitions(xk, uk);

    c(c_num*ct-c_num+1) = -xk1(1)+x_min;
    c(c_num*ct-c_num+2) = xk1(1)-x_max;
    c(c_num*ct-c_num+3) = -xk1(2)+x_min;
    c(c_num*ct-c_num+4) = xk1(2)-x_max;
    c(c_num*ct-c_num+5) = -uk(1)+u_min;
    c(c_num*ct-c_num+6) = uk(1)-u_max;
    c(c_num*ct-c_num+7) = -uk(2)+u_min;
    c(c_num*ct-c_num+8) = uk(2)-u_max;
 
    % update plant state and input for next step
    xk = xk1;
    if ct<N
        uk = u(:,ct+1);
    end
end
%% No equality constraints
ceq = [];

