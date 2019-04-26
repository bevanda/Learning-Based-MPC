function [cieq, ceq] = constraintsNMPC(dynamics,c,theta,x,N,xw,r0,K,state_F,state_h,in_F,in_h,term_F,term_h,Ts)
%% Constraint function for NMPC control of the Moore-Greitzer Compressor Model
%
% Inputs:
%   c, theta:      optimization variables
%   x:      current state at time k
%   N:      prediction horizon
%
% Outputs:
%   cineq:      inequality constraints applied across prediction horizon
%   ceq:    equality constraints 


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
        [xk1, uk] = dynamics(xk,ck,xw,r0,K,Ts);
        dxk=xk-xw;
        duk=uk-r0;
        dxk1=xk1-xw;
        % H-representation of constraints
        % state constraints
        cieq_run1 = state_F*dxk1-state_h;
        cieq  = [cieq; cieq_run1];
        % input constraints
        cieq_run2 = in_F*duk-in_h;
        cieq  = [cieq; cieq_run2];
        % update plant state and input for next step
        xk = xk1;
        if k<N
            ck = c(:,k+1);
        end
    else
        cieq_T = term_F*[dxk1;theta]-term_h;
%         cieq_T = [];
        cieq = [cieq; cieq_T];
    end
       
end

%% Equality constraints
ceq = [];

end
