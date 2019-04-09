function [xk, u] = nominalModel(x, u)
%% Discrete-time linear dynamic model 
%
% 4 states (x): 
%
% 1 input: (u)
%   
%
% 1 output: (y)
%   same as states (i.e. all the states are measureable)
%
% [A B C D] are state space matrices linearized at the current operating point.
%

%% Discrete time state-space model of the non-square LTI system for tracking
A = [1.01136382181963,0.0100343559666203,-6.46049734470989e-05,-1.93718915801510e-07; ...
    0.0100343559666203,0.995615459673586,-0.0127686112556342,-5.57236663816276e-05;...
    0,0,0.957038195891878,0.00792982548734094;...
    0,0,-7.92982548734093,0.602405619103784];
B = [-4.95341630475791e-07;-0.000193161656467182;0.0429618041081219;7.92982548734093];
C = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
% D = [0;0;0;0];
xk = A*x + B*u;
y = C*x;
end
