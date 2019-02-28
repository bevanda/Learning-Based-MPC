function [out1,out2, out3, out4] = solveMGCM(in)
%SOLVEMGCM Summary of this function goes here
%   Detailed explanation goes here
syms u ... % control input
    x1 ... % mass flow
    x2 ... % pressure rise
    x3 ... % throttle opeining
    x4 ... % throttle opening rate
    xk1 ...
    xk2 ...
    xk3 ...
    xk4 ...
    ;
wn=sqrt(1000); % resonant frequency
ksi=1/sqrt(2); % damping coefficient
beta=1; % constant >0
x2_c=0; % pressure constant

%% Continous time state-space model of the Moore-Greitzer compressor model

f1 = x2+x2_c+1+3*(x1/2)-(x1^3/2); % mass flow rate
f2 = (x1+1-x3*sqrt(x2))/beta; % pressure rise rate
f3 = x4; % throttle opening rate
f4 = -wn^2*x3-2*ksi*wn*x4+wn^2*in; % throttle opening acceleration

end

