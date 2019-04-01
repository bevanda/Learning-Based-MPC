function [xk1, uk] = getTransitions(xk, uk,sys)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
[xk1,~] = systemdynamics(xk, uk,sys);

end
