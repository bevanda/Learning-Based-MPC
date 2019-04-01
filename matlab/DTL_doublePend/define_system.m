function [A,B,C,Mtheta,Ntheta] = define_system()
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
D = [0 0];

n = size(A,1);
m = size(B,2);
o = size(C,1);
Ntheta = [1, 0];
% MN = [Mtheta; Ntheta];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
% Mtheta = [1, 0, 0, 0; 0, 1, 1, -2]';
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);
%% Initial disturbance estimate
d_0 = [0 0]';

end

