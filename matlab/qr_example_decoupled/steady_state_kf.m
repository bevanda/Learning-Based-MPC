function L = steady_state_kf(A, C, Q, R)

% A - system dynamics matrix, n x n
% C - output map, no x n
% Q - process covariance, nxn
% R - sensor covariance, no x no

[X,cl_eigs,G] = dare(A', C', Q, R);

L = G';