function [cieq, ceq] = constraintsFunction(c,theta,x,N,xw,r0,K,LAMBDA,PSI,state_F,state_h,in_F,in_h,term_F,term_h)
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
    if (k<N)
        % inequality constraints
        
        % obtain new state at next prediction step
        [xk1, uk] = getTransitionsTrue(xk, ck,xw,r0, K);
        x=xk1-xw;
        u=uk-r0;
%         Xs=E(1:n,:); Us= E(n+1:n+m);
        % H-representation of constraints
        % state constraints
        cieq_run1 = state_F*x-state_h;
        cieq  = [cieq; cieq_run1];
        % input constraints
        cieq_run2 = in_F*u-in_h;
        cieq  = [cieq; cieq_run2];
        % update plant state and input for next step
        xk = xk1;
        if k<N
            ck = c(:,k+1);
        end
    else
        cieq_T = term_F*[x;theta]-term_h;
%         cieq_T = [];
        cieq = [cieq; cieq_T];
    end
       
end

%% No equality constraints
ceq = [];

end
