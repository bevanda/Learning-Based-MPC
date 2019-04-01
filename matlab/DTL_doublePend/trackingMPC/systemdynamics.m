function [xk, y] = systemdynamics(x, u,sys)
%% Discrete time state-space model of the non-square LTI system for tracking
xk = sys.A*x + sys.B*u;
y = sys.C*x;
end
