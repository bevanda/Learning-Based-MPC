function [xk, y, A, B] = systemdynamics(x, u)
%% Continuous-time nonlinear dynamic model of a pendulum on a cart
%
% 2 states (x): 
%
% 2 inputs: (u)
%   
%
% 1 outputs: (y)
%   same as states (i.e. all the states are measureable)
%
% [A B C D] are state space matrices linearized at the current operating point.
%
% Copyright 2016 The MathWorks, Inc.

%#codegen
% Ntheta = [1, 0];
%% Discrete time state-space model of the non-square LTI system for tracking
A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1, 0];
xk = A*x + B*u;
y = C*x;
end
