function [cieq] = funLB(c,x,N,K,F,h)
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

%% Inequality constraints calculation
cieq = [];
% ceq = zeros(N*conseq_num,1);
% Apply cons_num*N state constraints across prediction horizon, from time
% k+1 to k+N
xk = x;
ck = c(:,1);
for k=1:N

    % inequality constraints

    % obtain new state at next prediction step
    [xk1, uk] = getTransitions(xk, ck, K);

    % H-representation of constraints
        cieq_temp = F*uk-h;
        cieq  = [cieq; cieq_temp];
    % update plant state and input for next step
    xk = xk1;
    if k<N
        ck = c(:,k+1);
    end

       
end


end
