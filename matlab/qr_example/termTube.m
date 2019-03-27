% Algorithm from:
% Theory and computation of disturbance invariant sets for discrete-time linear systems
% Ilya Kolmanovsky and Elmer G. Gilbert

%%%%%%%%%%%%%%%%% Calculation of the (R_N x {0}) set %%%%%%%%%%%%%%%%%%%%%%
function [F, h] = termTube(F_u, h_u, F_v, h_v)
% [F_u, h_u] -> MPI 
% [F_v, h_v] -> disturbance
% width of the terminal tube of the extended state for tracking: (R_N x {0}) 
N = length(h_u);
h = zeros(N,1);
options = optimset('Display', 'off');
% calculating the (R_N x {0}) set
for ind = 1:N
    % for each constraint calculate:
    % max_{x}(F_u(i,:)*x) s.t. F_v*x <= h_v
	[~, fval] = linprog(-F_u(ind,:),F_v,h_v,[], [], [], [], [], options);
	h(ind) = -fval;
end
F = F_u;
end