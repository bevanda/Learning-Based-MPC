% Algorithm from:
% Theory and computation of disturbance invariant sets for discrete-time linear systems
% Ilya Kolmanovsky and Elmer G. Gilbert

function [F_d, h_d] = pdiff(F_u, h_u, F_v, h_v)
options = optimset('Display', 'off');

N = length(h_u);
h = zeros(N,1);
for ind = 1:N
    % for each constraint calculate:
    %% max_{x}(F_u(i,:)*x) s.t. F_v*x <= h_v
	[~, fval] = linprog(-F_u(ind,:),F_v,h_v,[], [], [], [], [], options);
	h(ind) = -fval;
end
F_d = F_u;
h_d = h_u - h;