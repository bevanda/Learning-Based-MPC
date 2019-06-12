function cost = run_cost(x, u, theta, Q, R, LAMBDA, PSI)
    % Provide the running cost   
    cost = (x-LAMBDA*theta)'*Q*(x-LAMBDA*theta)+...
        (u-PSI*theta)'*R*(u-PSI*theta);
end

