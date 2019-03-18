function [cineq, ceq] = constraintsFunction(c,theta,x,N,K,LAMBDA)
%% Constraint function of nonlinear MPC for pendulum swing-up and balancing control
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%
% Output:
%   cineq:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints (empty)
%
% Copyright 2016 The MathWorks, Inc.

%% Nonlinear MPC design parameters
cons_num = 8;

x_min = -5.0; x_max=5.0;
u_min=-0.3; u_max=0.3;
%% Inequality constraints calculation
cineq = zeros(N*cons_num,1);
% Apply cons_num*N state constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
ck = c;

for k=1:N
    % obtain new state at next prediction step
    [xk1, uk] = getTransitions(xk, ck, K, theta, LAMBDA);

    cineq(cons_num*k-cons_num+1) = -xk1(1)+x_min;
    cineq(cons_num*k-cons_num+2) = xk1(1)-x_max;
    cineq(cons_num*k-cons_num+3) = -xk1(2)+x_min;
    cineq(cons_num*k-cons_num+4) = xk1(2)-x_max;
    cineq(cons_num*k-cons_num+5) = -uk(1)+u_min;
    cineq(cons_num*k-cons_num+6) = uk(1)-u_max;
    cineq(cons_num*k-cons_num+7) = -uk(2)+u_min;
    cineq(cons_num*k-cons_num+8) = uk(2)-u_max;
    
    % update plant state and input for next step
    xk = xk1;
    if k<N
        ck = c(:,k+1);
    end
       
end

%% No equality constraints
ceq = [];

end
