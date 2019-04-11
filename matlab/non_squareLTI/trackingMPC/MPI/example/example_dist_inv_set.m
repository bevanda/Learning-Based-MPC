addpath('../src/')
addpath('../src/utils/')

% make your own discrete linear system
A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
n = size(A,1);
m = size(B,2);
o = size(C,1);

Q = diag([1,1]);
R = diag([1,1]);
% K = -dlqr(A, B, Q, R);

mysys = LinearSystem(A, B, Q, R);

% constraints on state Xc and input Uc
% edges of the polytope
Xc_vertex = [2 -2;...
              2 2; ...
            -10 2;...
            -10 -2];
Uc_vertex = [1; -1];
F=[eye(2);-eye(2)];
h=[2;2;10;2];
Xc=polytope(F,h);
% Xc = Polyhedron(Xc_vertex);
figure
Xc.plot()
Uc = Polyhedron(Uc_vertex);
figure
Uc.plot()
%%
% construct a convex set of system noise (2dim here)
% W_vertex = [0.15, 0.15; 0.15, -0.15; -0.15, -0.15; -0.15, 0.15];
% W = Polyhedron(W_vertex);
% set boundary of system noise which corresponds to W.
w_max = [0.15; 0.15];
w_min = [-0.15; -0.15];
F_w=[eye(2);-eye(2)];
h_w=[w_max;-w_min];
W=Polyhedron(F_w,h_w);
figure
W.plot();
% compute disturvance invariant set Z.
Z = mysys.compute_distinv_set(W, 3, 1.05);

%%
% propagate particles many times following the discrete stochastic dynamics,
% and we will see that particles never go outside of distervance invariant set Z.
Nptcl = 10000;
x = zeros(2, Nptcl); % particles 
for i = 1:100
    sys_noise = rand(2, Nptcl).*repmat(w_max - w_min, 1, Nptcl) + repmat(w_min, 1, Nptcl);
    x = mysys.Ak*x + sys_noise;
    clf;
    Graphics.show_convex(Z, 'g', 'FaceAlpha', .3); % show Z
    scatter(x(1, :), x(2, :)); % show particles
    pause(0.01)
end
