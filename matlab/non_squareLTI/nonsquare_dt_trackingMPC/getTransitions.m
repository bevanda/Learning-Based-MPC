function [xk1, yk] = getTransitions(xk, uk)
%% Discrete-time nonlinear dynamic model of a pendulum on a cart at time k
%
% xk1 is the states at time k+1.
%

[xk1, yk] = systemdynamics(xk, uk);

end
