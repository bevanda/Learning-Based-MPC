function [xk1, uk, E] = getTransitions(xk, ck, K)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
r0=1.1547;
% K=-[+3.0741 2.0957 0.1197 -0.0090]; % K stabilising gain from the papers
K=[-4.53562790964753,-2.18725843827833,-0.116986431132625,0.00895660359632110];
uk = -K*xk+ck;
[xk1] = systemdynamics(xk, uk);

end
