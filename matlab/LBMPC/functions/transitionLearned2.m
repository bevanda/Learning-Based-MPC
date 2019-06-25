function [xk1, uk] = transitionLearned2(xk, ck, K, data)
%% Learned system transition after actuation
%
% Inputs:
%   xk: states at time k
%   ck: decision variable
%   K: stabilizing feedback gain
%   data: struct with X, Y obervations for estimation
%
% Outputs:
%   xk1: states at time k+1

uk = K*xk+ck;
[xk1] = learnedModel2(xk, uk, data);

end
