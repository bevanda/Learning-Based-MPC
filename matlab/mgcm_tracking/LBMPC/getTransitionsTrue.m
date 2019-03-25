function [xk1, uk, E] = getTransitionsTrue(xk, ck, K)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
r0=1.1547;
% K=[+3.0741 2.0957 0.1197 -0.0090]; % K stabilising gain from the papers
K=[3.0742   -2.0958   -0.1194    0.0089];
uk = -K*xk+r0;
[xk1] = trueModel(xk, uk);

end
