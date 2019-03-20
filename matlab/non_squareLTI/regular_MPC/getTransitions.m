function [xk1, uk] = getTransitions(xk, c)
%% Discrete-time dynamic plant with prestabilisation
%
% xk1 is the states at time k+1.
% c - decision variable
uk = c;
[xk1,~,~,~] = systemdynamics(xk, uk);

end
