function [xk1, uk, E] = getTransitions(xk, ck, K)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.

uk =K*xk+ck;
[xk1] = systemdynamics(xk, uk);

end
