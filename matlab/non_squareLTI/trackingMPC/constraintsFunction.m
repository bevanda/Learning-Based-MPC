function [cieq, ceq] = constraintsFunction(c,theta,x,N,K,LAMBDA,PSI)
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
consieq_num = 16;
% conseq_num  = 8;
ALPHA = 0.99; % λ ∈ (0, 1), λ can be chosen arbitrarily close to 1, the obtained
% invariant set can be used as a reliable polyhedral approximation to the maximal invariant set 

x_min = -5.0; x_max=5.0;
u_min=-0.3; u_max=0.3;
%% Inequality constraints calculation
cieq = zeros(N*consieq_num,1);
% ceq = zeros(N*conseq_num,1);
% Apply cons_num*N state constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
ck = c;
n=size(xk,1); m=size(ck,1);
for k=1:N
    % obtain new state at next prediction step
    [xk1, uk, E] = getTransitions(xk, ck, K, theta, LAMBDA,PSI);
    Xs=E(1:n,:); Us= E(n+1:n+m);
    % inequality constraints
    cieq(consieq_num*k-consieq_num+1) = -xk1(1)+x_min;
    cieq(consieq_num*k-consieq_num+2) = xk1(1)-x_max;
    cieq(consieq_num*k-consieq_num+3) = -xk1(2)+x_min;
    cieq(consieq_num*k-consieq_num+4) = xk1(2)-x_max;
    cieq(consieq_num*k-consieq_num+5) = -uk(1)+u_min;
    cieq(consieq_num*k-consieq_num+6) = uk(1)-u_max;
    cieq(consieq_num*k-consieq_num+7) = -uk(2)+u_min;
    cieq(consieq_num*k-consieq_num+8) = uk(2)-u_max;
    if k==N
    % terminal equality constraints
        cieq(consieq_num*k-consieq_num+9) = -Xs(1)+x_min*ALPHA;
        cieq(consieq_num*k-consieq_num+10) = Xs(1)-x_max*ALPHA;
        cieq(consieq_num*k-consieq_num+11) = -Xs(2)+x_min*ALPHA;
        cieq(consieq_num*k-consieq_num+12) = Xs(2)-x_max*ALPHA;
        cieq(consieq_num*k-consieq_num+13) = -Us(1)+u_min*ALPHA;
        cieq(consieq_num*k-consieq_num+14) = Us(1)-u_max*ALPHA;
        cieq(consieq_num*k-consieq_num+15) = -Us(2)+u_min*ALPHA;
        cieq(consieq_num*k-consieq_num+16) = Us(2)-u_max*ALPHA;
    end
    % update plant state and input for next step
    xk = xk1;
    if k<N
        ck = c(:,k+1);
    end
       
end

%% No equality constraints
ceq = [];

end
