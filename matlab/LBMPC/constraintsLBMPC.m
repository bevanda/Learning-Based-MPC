function [cieq, ceq] = constraintsLBMPC(c,theta,x,N,K,state_F,state_h,in_F,in_h,term_F,term_h,stateD_F,stateD_h)
%% Constraint function for LBMPC control of the Moore-Greitzer Compressor Model
%
% Inputs:
%   c, theta:      optimization variables
%   x:      current state at time k
%   N:      prediction horizon
%
% Outputs:
%   cineq:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints 
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
        % obtain new state at next prediction step
        [xk1, uk] = transitionNominal(xk, ck, K);
        % H-representation of constraints
        % contracted state constraint for every the k+1 step
        if k == 1
            cieq_run = stateD_F*xk1-stateD_h;
            cieq  = [cieq; cieq_run]; %#ok<*AGROW>
        end
        % state constraints
        cieq_run1 = state_F*xk1-state_h;
        cieq  = [cieq; cieq_run1]; %#ok<*AGROW>
        % input constraints
        cieq_run2 = in_F*uk-in_h;
        cieq  = [cieq; cieq_run2]; %#ok<*AGROW>
        % update plant state and input for next step
        xk = xk1;
        if k<N
            ck = c(:,k+1);
        end
    else
        cieq_T = term_F*[xk1;theta]-term_h;
        cieq = [cieq; cieq_T]; %#ok<*AGROW>
    end
       
end
%% Equality constraints
ceq = [];

end
