function [A,B,C,d_0] = define_system()
%% Simple discrete non-square LTI system 

%==========================================================================
%% Define A, B, C of discrete time model:
% x+ = Ax + Bu
% y = Cx
%
%==========================================================================
A = [1 ,1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];

%% Initial disturbance estimate
d_0 = [0 0]';

end

