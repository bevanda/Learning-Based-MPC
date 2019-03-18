function [xk1, uk] = getTransitions(xk, c, K, theta, LAMBDA)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
% c - decision variable
uk = -K*xk+LAMBDA*theta+c;
[xk1] = systemdynamics(xk, uk);

end
