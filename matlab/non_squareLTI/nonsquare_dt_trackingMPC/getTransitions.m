function [xk1, yk] = getTransitions(xk, uk, K, xs)
%% Discrete-time dynamic model at time k
%
% xk1 is the states at time k+1.
%

[xk1, yk] = systemdynamics(xk, -K*(xk-xs)+uk);

end
