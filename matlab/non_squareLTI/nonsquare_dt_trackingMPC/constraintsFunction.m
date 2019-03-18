function [c, ceq] = constraintsFunction(u,theta,x,N,K,xs)
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
cons_num = 4;

x_min = -inf; x_max=5.0;
u_min=-inf; u_max=0.3;
%% Inequality constraints calculation
c = zeros(N*cons_num*2,1);
% Apply cons_num*N state constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
uk = u;

for k=1:N
    % obtain new state at next prediction step
    xk1 = getTransitions(xk, uk, K, xs);

%     c(cons_num*k-cons_num+1) = -xk1(1)+x_min;
    c(cons_num*k-cons_num+2) = xk1(1)-x_max;
%     c(cons_num*k-cons_num+3) = -xk1(2)+x_min;
    c(cons_num*k-cons_num+4) = xk1(2)-x_max;
%     c(cons_num*k-cons_num+5) = -uk(1)+u_min;
    c(cons_num*k-cons_num+6) = uk(1)-u_max;
%     c(cons_num*k-cons_num+7) = -uk(2)+u_min;
    c(cons_num*k-cons_num+8) = uk(2)-u_max;
%     c(cons_num*k-cons_num+9) = -theta(1)+inf;
%     c(cons_num*k-cons_num+10) = theta(1)-inf;
%     c(cons_num*k-cons_num+11) = -theta(2)+inf;
%     c(cons_num*k-cons_num+12) = theta(2)-inf;
    % update plant state and input for next step
    xk = xk1;
    if k<N
        uk = u(:,k+1);
    end
       
end
% c_theta = zeros(4,1);
% if k==N
%     c_theta(1) = -theta(1)+inf;
%     c_theta(2) = theta(1)-inf;
%     c_theta(3) = -theta(2)+inf;
%     c_theta(4) = theta(2)-inf;
%     c = [c; c_theta];
% end

%% No equality constraints
ceq = [];

end
