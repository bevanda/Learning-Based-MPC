function [xk1, uk, E] = getTransitions(xk, ck,xw,r0, K)
%% Discrete-time dynamic plant with prestabilisation
%

uk = K*(xk-xw)+ck+r0;
[xk1] = systemdynamics(xk, uk);

end
