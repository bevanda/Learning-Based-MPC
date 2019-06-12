function cost = cost_fun(N, y, X, x_eq, A,B,Q, R, P, T, m, LAMBDA, PSI)
    % Formulate the cost function to be minimized
    cost = 0;
    u=y(1:m*N);
    theta=y(m*N+1:end);
    % k+1 to k+N
    xk = X;
    uk = u(1:m);
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N
        % obtain new state at next prediction step
        xk1 = dynamics(xk, uk, A, B);
        cost = cost + run_cost(xk, uk, theta, Q, R, LAMBDA,PSI);
        % Update xk and uk for the next prediction step.
        xk = xk1;
        uk = u(m*k+1:m*(k+1));
    end
    cost = cost + term_cost(x, theta, x_eq, P, T, LAMBDA);
end