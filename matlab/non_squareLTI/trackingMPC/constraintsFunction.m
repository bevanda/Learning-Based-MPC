function [cieq, ceq] = constraintsFunction(u,theta,x,N,K,LAMBDA,PSI,run_F, run_h, term_F, term_h)
%% Constraint function of MPC 
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

%% Inequality constraints calculation
cieq = [];
% ceq = zeros(N*conseq_num,1);
% Apply cons_num*N state constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
uk = u(:,1);
n=size(xk,1); m=size(uk,1);
for k=1:N

    % obtain new state at next prediction step
    [xk1, uk] = getTransitions(xk, uk);
    cieq = [cieq; run_F*[xk;uk]-run_h];
    % update plant state and input for next step
    xk = xk1;
    if k<N
        uk = u(:,k+1);
    end
    if k==N
        cieq_T = term_F*[xk;theta]-term_h;
        cieq = [cieq; cieq_T];
    end
       
end

%% No equality constraints
ceq = [];

end
