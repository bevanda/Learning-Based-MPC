function J = NMPC_costFunction(c,theta,x,xs,N,c0,Q,R,P,T,K,x_w,r0,LAMBDA,PSI,data,dT)
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
%

%% LMPC design parameters

% Set initial plant states, controller output and cost.
xk = x;
ck = c0;

J = 0;
% Loop through each prediction step.
for k=1:N-1  
    % Obtain plant state at next prediction step.
%     [xk1,uk]= getTransitionsLearn(xk,ck,K,data); % learned model
%     [xk1,uk]=getTransitions(xk,ck,K);
    [xk1, uk] = getTransitionsTrue(xk,ck,x_w,r0,K,dT); % plant  
    dx=xk-x_w;
    du=uk-r0;
    % RUNNING COST
    % accumulate state tracking cost from x(k+1) to x(k+N).
    J = J + (dx-LAMBDA*theta)'*Q*(dx-LAMBDA*theta);
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    J = J + (du-PSI*theta)'*R*(du-PSI*theta);

    % Update xk and uk for the next prediction step.
    xk = xk1;
 
    ck = c(:,k+1);

end
%TERMINAL COST
J = J + (dx-LAMBDA*theta)'*P*(dx-LAMBDA*theta) + (LAMBDA*theta-xs)'*T*(LAMBDA*theta-xs);
end
