function [xk1, uk, E] = getTransitionsLearn(xk, ck, K,data)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.

K=[3.0742   -2.0958   -0.1194    0.0089];
uk = -K*xk+ck;
[xk1] = learnedDynamics(xk, uk,data);

end
