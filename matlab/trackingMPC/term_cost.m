function [cost] = term_cost(x, theta, x_eq, P, T,LAMBDA)
    % Terminal cost
    cost = (x-LAMBDA*theta)'*P*(x-LAMBDA*theta)+(LAMBDA*theta - x_eq)'*T*(LAMBDA*theta - x_eq);
end