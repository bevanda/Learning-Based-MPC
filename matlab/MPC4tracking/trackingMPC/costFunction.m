function J = costFunction(c,theta,x,xs,N,u0,P,T,K,LAMBDA,PSI,sys)
%% Cost function of non-quadratic discrete time LTI
% Inputs:
%   c:      decision variable, from time k to time k+N-1 
%   theta:  decision variable, state and input parametrisation 
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%   xref:   state references, constant from time k+1 to k+N
%   u0:     previous controller output at time k-1
%   theta0: previous theta output at time k-1

% Output:
%   J:      objective function cost

Q = sys.Q;
R = sys.R;

% Set initial plant states, controller output and cost
xk = x;
uk = u0;
J = 0;
% Loop through each prediction step.
%%%%% RUNNING COST %%%%%
for k=1:N-1  
    % accumulate state tracking cost from x(k+1) to x(k+N).
    J = J + (xk-LAMBDA*theta)'*Q*(xk-LAMBDA*theta);
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    J = J + (uk-PSI*theta)'*R*(uk-PSI*theta);
    
    % Obtain plant state at next prediction step.
    [xk1,~]= getTransitions(xk, uk, sys);
    % Update xk and uk for the next prediction step.
    xk = xk1;
    uk = c(:,k+1);
end
%%%%% TERMINAL COST %%%%%
J = J + (xk-LAMBDA*theta)'*P*(xk-LAMBDA*theta) + ...
    (LAMBDA*theta-xs)'*T*(LAMBDA*theta-xs);
end


