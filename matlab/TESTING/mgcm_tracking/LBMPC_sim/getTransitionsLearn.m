function [xk1, uk] = getTransitionsLearn(xk, ck, K,data)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
uk = K*xk+ck;
[xk1] = learnedDynamics(xk, uk, data);

end
