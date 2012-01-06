%==========================================================================
% Configuration for MPT toolbox
%==========================================================================
if strcmp(getenv('USER'), 'bouffard') == 0
    % Stuff for Anil only
    mpt_path = 'C:\Program Files\MATLAB\R2009b\toolbox\mpt';
    disp('Adding MPT toolbox to MATLAB path...');
    addpath(genpath(mpt_path));
    
    %% Misc configuration
    quad_dat_fname = 'quad.dat';
    quad_dat_fname = 'quad.mat';
    dual_ekf_fname = 'dual_ekf.dat';
    
    % added by George
    % quad_bin_fname = 'quad.bin';
else
    % Stuff for Pat only
    clear all;
    %clc;
    close all;
    
    %% Misc configuration
    quad_bin_fname = 'quad.bin';
end

%% Get system definition:
add_alt_dynamics = 1;
[A,B,C,d_0] = define_system(add_alt_dynamics);
% Count number of states n, number of inputs m, number of outputs o:
n = size(A,1);
m = size(B,2);
o = size(C,1);
q = 5; % number of outputs (for dual EKF purposes .. should o also be set to 4?)


%% Get design parameters:
[N, Q, R, ...
    dlqr_controlweight, maxadm_controlweight, ...
    max_x, max_y, max_vx, max_vy, ...
    max_z, min_z, max_vz, ...
    max_pitch_cmd, max_roll_cmd, ...
    max_thrust_cmd, min_thrust_cmd, ...
    max_pitch, max_pitch_rate, max_roll, max_roll_rate, ...
    state_uncert, ...
    enable_learning, ALPHA, MU_factor, ...
    uncertainty_block] = design_params(n, m);
%==========================================================================

%%
% feedback policy K_t used for terminal set computations:
K_t = -dlqr(A, B, Q, maxadm_controlweight*R);

%==========================================================================
% Define a nominal feedback policy K and corresponding terminal cost
% Feedback policy chosen using discrete time linear quadratic regulator.
% This can be changed based on engineering preference:
K = -dlqr(A, B, Q, dlqr_controlweight*R);

% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, dlqr_controlweight*R);

%%
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================

temp0 = [max_pitch_cmd; max_roll_cmd; max_thrust_cmd];
temp01 = [max_pitch_cmd; max_roll_cmd; -min_thrust_cmd];
temp0 = temp0(1:m);
temp01 = temp01(1:m);
F_u = [eye(m); -eye(m)]; h_u = [temp0;temp01];
temp1 = [max_x; max_vx; max_pitch; max_pitch_rate; max_y; max_vy; max_roll; max_roll_rate; max_z; max_vz];
temp2 = [max_x; max_vx; max_pitch; max_pitch_rate; max_y; max_vy; max_roll; max_roll_rate; -min_z; max_vz];
F_x = [eye(n); -eye(n)]; h_x = [temp1;temp2];
F_g = [eye(n); -eye(n)]; h_g = [state_uncert; state_uncert]; % uncertainty polytope

% count the length of the constraints on input, states, and uncertainty:

length_Fu = length(h_u);
length_Fx = length(h_x);
length_Fg = length(h_g);

%% State constraints:
% Points that after (A+B*K_t) get to (F_x,h_x) \ominus (F_g,h_g)
[F_x_g, h_x_g] = double(polytope(F_x, h_x) - polytope(F_g, h_g));
Fx{1} = F_x_g;
fx{1} = h_x_g;
for i=2:N
    Fx{i} = F_x;
    fx{i} = h_x;
end
for i=1:N
   Fu{i} = F_u;
   fu{i} = h_u;
end




%%
%==========================================================================
% Generate steady-state parameterization and their projection matrices
%==========================================================================

M = [A - eye(n), B, zeros(n,o); ...
    C, zeros(o,m), -eye(o)];
V = null(M);
LAMBDA = V(1:n,:);
PSI = V(n+1:n+m,:);
XI = V(n+m+1:n+m+o,:);

% Solutions of M*[x;u;y] = [-d;0] are of the form
% M\[-d;0] + V*theta, theta in R^m
V_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = V_0(1:n);
PSI_0 = V_0(n+1:n+m);
XI_0 = V_0(n+m+1:n+m+o);

% these 2 are not used but we'll keep 'em for now..
Proj_X_t = LAMBDA*pinv(LAMBDA);
Proj_U_t = PSI*pinv(LAMBDA);

%%
%==========================================================================
% Compute maximally invariant set
%==========================================================================

disp('Computing and simplifying terminal set...');
F_w = [F_x zeros(length_Fx, m);
    zeros(length_Fx, n) F_x*LAMBDA; ...
    F_u*K_t, F_u*(PSI - K_t*LAMBDA); ...
    zeros(length_Fu, n) F_u*PSI; ...
    F_x_g*(A+B*K_t) F_x_g*B*(PSI-K_t*LAMBDA)];
h_w = [...
    h_x; ...
    h_x - F_x*LAMBDA_0; ...
    h_u - F_u*(PSI_0 - K_t*LAMBDA_0); ...
    h_u - F_u*PSI_0; ...
    h_x_g - F_x_g*B*(PSI_0-K_t*LAMBDA_0)];




% Subtract out points due to disturbance (F_g,h_g)
[F_w_N0, h_w_N0] = pdiff(F_w, h_w, ...
    [F_g zeros(length_Fg,m); ...
    zeros(m, n) eye(m); ...
    zeros(m, n) -eye(m)], ...
    [h_g; ...
    zeros(2*m,1)]);

% Simplify the constraints
term_poly = polytope(F_w_N0, h_w_N0);
[F_w_N, h_w_N] = double(term_poly);

disp('Terminal set polytope:');
term_poly

% x-theta constraints:
F_xTheta = F_w(:, 1:n);
F_theta = F_w(:,n+1:n+m);
f_xTheta = h_w;

%%
%==========================================================================
% Generate inequality constraints
%==========================================================================

length_Fw = size(F_w_N, 1);

Aineq = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, N*m+m);
bineq = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, 1);
b_crx = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, n);

L_i = zeros(n, N*m);
KL_i = zeros(m, N*m);
disp('Generating constraints on inputs...');
d_i = zeros(n,1);
for ind = 1:N
    %disp(['u ind: ', num2str(ind)]);
    
    KL_i = K*L_i;
    KL_i(:, (ind-1)*m + (1:m)) = eye(m);
    
    Aineq((ind-1)*length_Fu + (1:length_Fu),1:N*m) = F_u*KL_i;
    bineq((ind-1)*length_Fu + (1:length_Fu)) = h_u - F_u*K*d_i;
    b_crx((ind-1)*length_Fu + (1:length_Fu),:) = -F_u*K*(A+B*K)^(ind-1);
    
    L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];
    d_i = (A+B*K)*d_i + d_0;
end

L_i = zeros(n, N*m);
disp('Generating constraints on states...');
d_i = d_0;
for ind = 1:N
    %disp(['x ind: ', num2str(ind)]);
    L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];
    
    if ind == 1
        disp('Generating terminal constraints on states...');
        Aineq(N*length_Fu + (1:length_Fw), :) = F_w_N*[L_i zeros(n,m); zeros(m,N*m) eye(m)];
        bineq(N*length_Fu + (1:length_Fw)) = h_w_N - F_w_N*[d_i; zeros(m,1)];
        b_crx(N*length_Fu + (1:length_Fw),:) = -F_w_N*[(A+B*K)^(ind); zeros(m,n)];
        
    else
        
        Aineq(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx),1:N*m) = F_x*L_i;
        bineq(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx)) = h_x - F_x*d_i;
        b_crx(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx),:) = -F_x*(A+B*K)^(ind);
        
    end
    
    d_i = (A+B*K)*d_i + d_0;
end

ind = N;
L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];


CONSTRAINT_COUNT = length(bineq);
% disp('Removing redundant inequality constraints...');
% temp_tope = polytope(Aineq, bineq);
% [Aineq, bineq] = double(temp_tope);

%%
%==========================================================================
% Generate uncertainty structure format for the C++ code
%==========================================================================

disp('Generating uncertainty structure...');
sub_block = sum(uncertainty_block, 2);
uncertainty_structure = zeros(n, n + m + 1);
for ind = 1:n
    uncertainty_structure(ind, 1:sub_block(ind)) = find(uncertainty_block(ind,:) == 1) - 1;
end

%%
%==========================================================================
% Compute uncertainty bounds
%==========================================================================

disp('Generating uncertainty bounds...');
options = optimset('Display', 'off');
uncertainty_upper_bounds = zeros(n, n + m + 1);
uncertainty_lower_bounds = zeros(n, n + m + 1);

X_vertices = extreme(polytope(blkdiag(F_x, F_u), [h_x; h_u]));
length_X_vertices = size(X_vertices,1);
X_vertices = [X_vertices ones(length_X_vertices,1)];

for ind = 1:n
    disp(['sub_block ind: ', num2str(ind)]);
    
    if (sub_block(ind) > 0)
        s_vec = zeros(n,1); s_vec(ind) = 1;
        F_g_s = [1; -1]; h_g_s = [s_vec'*linprog(-s_vec, F_g, h_g); -s_vec'*linprog(s_vec, F_g, h_g, [], [], [], [], [], options)];
        F_delta_A = kron(F_g_s, X_vertices(:, logical(uncertainty_block(ind, :))));
        h_delta_A = repmat(h_g_s, length_X_vertices, 1);
        [F_delta_A, h_delta_A] = double(polytope(F_delta_A, h_delta_A));
        
        for ind_j = 1:sub_block(ind)
            s_vec = zeros(sub_block(ind),1); s_vec(ind_j) = 1;
            uncertainty_upper_bounds(ind, uncertainty_structure(ind, ind_j) + 1) = s_vec'*linprog(-s_vec, F_delta_A, h_delta_A, [], [], [], [], [], options);
            uncertainty_lower_bounds(ind, uncertainty_structure(ind, ind_j) + 1) = s_vec'*linprog(s_vec, F_delta_A, h_delta_A, [], [], [], [], [], options);
        end
    else
        uncertainty_lower_bounds(ind, :) = 0;
        uncertainty_upper_bounds(ind, :) = 0;
    end
end

%%
% MU_VEC becomes the diagonal entries of MU, this is the weighting on the
% priors
if enable_learning
    MU_VEC = [MU_factor*ones(1,n+m+1)]; % weight on prior model (0 = completely uncertain)
    MU_VEC(13) = MU_VEC(13)*0.01;
else
    MU_VEC = [1e10*ones(1,n+m+1)]; % weight on prior model -- 'turn off' learning
end

%%
% Dual EKF stuff
p = 12; % number of parameters being estimated
P2_0 = zeros(n,p); % initial cross-covariance
P3_0 = zeros(p,p); % initial parameter covariance
%XI = diag([1e-4 4e-4 1e-4 4e-4 1e-4])
%XI = 0.01*eye(q); % innovation noise covariance

Cekf = zeros(q, n);
% we observe x pos, pitch, y pos, roll, and z pos:
Cekf(1,1) = 1;
Cekf(2,3) = 1;
Cekf(3,5) = 1;
Cekf(4,7) = 1;
Cekf(5,9) = 1;
% EKF feedback gain... how to pick??
sys = ss(A, B, Cekf, 0, -1); % -1 indicates discrete time plant w/ unspec. sample time

%[kest,Kekf,Pkf] = kalman(sys,1*eye(m),XI);
%Kekf = steady_state_kf(A, Cekf, eye(n), eye(q));
%Kekf = place(A', Cekf', [0.25 0.25 0.25 0.25 0.25 0.4 0.4 0.4 0.4 0.4])';
%Kekf = place(A', Cekf', [0.5 0.5 0.5 0.5 0.5 0.8 0.8 0.8 0.8 0.8])';
%Kekf = zeros(n,q);
%Kekf = steady_state_kf(A, Cekf, 1e-2*eye(n), XI); % KF too gentle --> params too aggressive
%Kekf = steady_state_kf(A, Cekf, 1e2*eye(n), 10000*XI); % KF too gentle --> params too aggressive


%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.6 0.8])';
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.8 0.9])'; % too fast - parameters unstable
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.7 0.8])'; % still too slow
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.75 0.85])'; % too fast - parameters unstable

% Unstable:
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.7 0.8])'; %
%XI = diag([0.01 0.01 0.01 0.01 0.001]);

% Stable:
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.8 + 0.1i 0.8 - 0.1i])'; %
%XI = diag([0.01 0.01 0.01 0.01 0.01]);

% Stable:
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.8 + 0.4i 0.8 - 0.4i])'; %
%XI = diag([0.01 0.01 0.01 0.01 0.01]);

% Unstable:
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.9 + 0.1i 0.9 - 0.1i])'; %
%XI = diag([0.01 0.01 0.01 0.01 0.01]);

% Stable:
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.85 + 0.4i 0.85 - 0.4i])'; %
%XI = diag([0.01 0.01 0.01 0.01 0.01]);

% Unstable:
%Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.875 + 0.5i 0.875 - 0.5i])'; %
%XI = diag([0.01 0.01 0.01 0.01 0.01]);

% Stable:
Kekf_x = place(A(1:4,1:4)', Cekf(1:2,1:4)', [0.25 0.25 0.4 0.4])';
Kekf_y = Kekf_x;
Kekf_z = place(A(9:10,9:10)', Cekf(5, 9:10)', [0.85 + 0.4i 0.85 - 0.4i])'; %
XI = diag([0.01 0.01 0.01 0.01 0.0075]);


Kekf = blkdiag(Kekf_x,Kekf_y,Kekf_z);

% Parameter bounds, ummm...
%foo = 0.2; % up to 20% error on A and B entries
%beta_min = zeros(p, 1);
nom_AB = [...
    0.05*A(4,3) 0.05*A(4,4) 0.05*B(4,1) ... % x-axis
    0.05*A(8,7) 0.05*A(8,8) 0.05*B(8,2) ... % y-axis
    0.15*B(10,3) ... % z-axis (ground effect)
    ]; %
%nom_AB = [A(4,3) A(4,4) B(4,1) A(8,7) A(8,8) B(8,2) 0.1*B(10,3)]; % ... but 5% on vert. B entries
beta_min = -abs([nom_AB zeros(1,5)])';
dT = 0.025; % timestep
max_lat_accel_noise = 1.0; % m/s^2
max_vert_accel_noise = 1.0*max_lat_accel_noise;
max_ang_accel_noise = 1.0; %rad/s^2
beta_min(8:end) = -[dT*max_lat_accel_noise dT*max_ang_accel_noise dT*max_lat_accel_noise dT*max_ang_accel_noise 0.5*dT^2*max_vert_accel_noise];
%beta_max = -beta_min;
%beta_min = -0.5*max_vert_accel_noise*dT^2;
beta_max = - beta_min;
P3_0 = diag(0.001*[...
    2 4 2 ... % x (b1, b2, b3)
    2 4 2 ... % y (b4, b5, b6)
    10 ... % z (b7)
    1e-3 ... % x pos b8
    1e-3 ... % x angle b9
    1e-3 ... % y pos b10
    1e-3 ... % y angle b11
    1e-3 ... % z pos b12
    ]'.*beta_max);
P3_lambda = 0.1*P3_0; %1e-4*diag(beta_max);

beta_max(7) = beta_max(7)*0.1; % only allow 1/10 the correction to this parameter in the positive sense

% ground_effect_fudge1 = 10; % nominal
% %ground_effect_fudge1 = 100; % 'bad learning'
% ground_effect_fudge2 = 0.001;
% P3_0(7,7) = ground_effect_fudge1 * P3_0(7,7);
% P3_lambda(7,7) = ground_effect_fudge1 * P3_lambda(7,7);
% P3_0(12,12) = ground_effect_fudge2 * P3_0(12,12);
% P3_lambda(12,12) = ground_effect_fudge2 * P3_lambda(12,12);

% 'bad learning'
%P3_0(2,2) = 5*P3_0(2,2);
%P3_lambda(2,2) = 5*P3_lambda(2,2);

%% Parameters for constructor

n_iter = 200; % maximum number of Newton iterations
reg = 1e-2;  % regularization Term
eps_primal = 0.1;
eps_dual = 0.1;
eps_mu = 0.1;

%% Write quad.dat and quad.mat
s = d_0;
Q_tilde = Q;
Q_tilde_f = Q;
writeParam(...
    N, ...
    quad_bin_fname, ...
    n_iter, ...
    reg, ...
    eps_primal, ...
    eps_dual, ...
    eps_mu, ...
    A, B, s, ...
    Q_tilde, ...
    Q_tilde_f, ...
    R, ...
    Fx, ...
    fx, ...
    Fu, ...
    fu, ...
    F_xTheta, ...
    F_theta, ...
    f_xTheta, ...
    K ...
    )

