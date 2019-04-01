function [xk1, uk] = getTransitions(xk, uk)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
[xk1,~,~,~] = systemdynamics(xk, uk);

end
