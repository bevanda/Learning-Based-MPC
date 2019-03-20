function [xk1, uk, E] = getTransitions(xk, ck, K, theta, LAMBDA, PSI)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
% c - decision variable
% K=[-3.0741 2.0957 0.1197 -0.0090]; % K stabilising gain from the papers
K= [-3.69438241357512	-1.78296525690465	-0.271009055947981	0.00595919621278219];
s=size(K,1);
Kw = [K eye(s)];
E =[LAMBDA' PSI']'*theta;
% r0=1.1547;
% uk = -K*xk+Kw*E+ck;
% uk = -K*xk+ck;
uk = ck;
[xk1] = systemdynamics(xk, uk);

end
