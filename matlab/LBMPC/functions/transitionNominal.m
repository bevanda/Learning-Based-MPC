function [xk1, uk] = transitionNominal(xk, ck, K)
%% Nominal system transition after actuation
%
% Inputs:
%   xk: states at time k
%   ck: decision variable
%   K: stabilizing feedback gain
%
% Outputs:
%   xk1: states at time k+1
%
uk = K*xk+ck;
[xk1] = nominalModel(xk, uk);

end
