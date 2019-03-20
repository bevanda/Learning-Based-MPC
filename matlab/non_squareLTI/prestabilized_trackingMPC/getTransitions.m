function [xk1, uk, E] = getTransitions(xk, ck, K, theta, LAMBDA, PSI)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
% c - decision variable
% size(xk,2)
m=size(K,2);
Kw = [-K eye(m)];
E =[LAMBDA' PSI']'*theta;
uk = K*xk+Kw*E+ck;
[xk1] = systemdynamics(xk, uk);
% size(xk1,2)
end
