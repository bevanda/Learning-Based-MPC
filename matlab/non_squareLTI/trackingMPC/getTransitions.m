function [xk1, uk, E] = getTransitions(xk, c, K, theta, LAMBDA, PSI)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
% c - decision variable
E =[LAMBDA' PSI']'*theta;
uk = c;
[xk1] = systemdynamics(xk, uk);

end
