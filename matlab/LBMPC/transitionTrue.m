function [xk1, uk] = transitionTrue(xk,ck,xw,r0,K,dT)
%% True system transition after actuation
% % Inputs:
%   xk: states at time k
%   ck: decision variable
%   K: stabilizing feedback gain
%
% Outputs:
%   xk1: states at time k+1
%
uk = K*(xk-xw)+ck+r0;
[xk1] = trueModel(xk, uk, dT);

end
