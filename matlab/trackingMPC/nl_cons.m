function [cieq] = nl_cons(N,y,X,m,A,B,F_run,h_run, F_term, h_term) 
% Introduce the nonlinear constraints also for the terminal state
   
u=y(1:m*(N));
theta=y(m*(N+1):end);

cieq = [];
% ceq = zeros(N*conseq_num,1);
% Apply cons_num*N state constraints across prediction horizon, from time
% k+1 to k+N
xk = X;
uk = u(1:m);
for k=1:N
    % obtain new state at next prediction step
    xk1 = dynamics(xk, uk,A,B);
    cieq = [cieq; F_run*[xk;uk]-h_run];
    % update plant state and input for next step
    xk = xk1;
    if k<N
        uk = u(m*k+1:m*(k+1));
    end
    if k==N
        cieq_T = F_term*[xk;theta]-h_term;
        cieq = [cieq; cieq_T];
    end
end

