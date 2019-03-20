function [xk, y, A, B] = systemdynamics(x, u)
%% Discrete-time linear dynamic model 
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
A = [1.01125000000000,0.0100000000000000,0,0;...
    0.0100000000000000,0.995555557627778,-0.0129903810567666,0;...
    0,0,1,0.0100000000000000;...
    0,0,-10,0.552786404500042];
B = [0;0;0;10];
C = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
% D = [0;0;0;0];
xk = A*x + B*u;
y = C*x;
end
