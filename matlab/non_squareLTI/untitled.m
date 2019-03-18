%% Non-square LTI system 
% Parameters
% Set the prediction horizon:
N = 3;

%% Closed-Loop Simulation
% The initial conditions
x = [0;...
    -2];
%setpoint
xs = [4.95;...
      0.0];

u_min =-0.3; u_max= 0.3;

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];
%% Configuration
mgcm_bin_fname = 'nonsq.bin';

N_values = [3 10 ... % horizons for which to calculate the discriminating kernel (maximum disturbance invarinat set)
    % 60 120 240 ...
    ];

for N=N_values
    nonsq_bin_fname = ['nonsq_N' num2str(N) '.bin'];
    %% Get system definition:
    [A,B,C,d_0] = define_system();
    %% Get design parameters:
    [Q, R, ...
    x_max, x_min, ...
    u_max, u_min, ...
    state_uncert, ...
    ALPHA, MU_factor] = design_params(n, m);
    %==========================================================================
    %%
    % Count number of states n, number of inputs m, number of outputs o:
    n = size(A,1);
    m = size(B,2);
    o = size(C,1);

    [P, e, K] = dare(A,B,Q,R);

    %%
    %==========================================================================
    % Define polytopic constraints on input F_u*x <= h_u and
    % state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
    %==========================================================================

    F_u = [eye(m); -eye(m)]; h_u = [u_max;u_max; -u_min;-u_min];
    F_x = [eye(n); -eye(n)]; h_x = [x_max;x_max; -x_min;-x_min];

    % count the length of the constraints on input, states, and uncertainty:
    length_Fu = length(h_u);
    length_Fx = length(h_x);


    %% State constraints:
    % Points that after (A+B*K) get to (F_x,h_x) \ominus (F_g,h_g)
%     [F_x_g, h_x_g] = double(polytope(F_x, h_x) - polytope(F_g, h_g));
%     tempPoly = Polyhedron(F_x, h_x) - Polyhedron(F_g, h_g);
%     F_x_g = tempPoly.A; % Inequality description { x | H*[x; -1] <= 0 }
%     h_x_g = tempPoly.b; % Inequality description { x | A*x <= b }
    
    for i=1:N
        Fx{i} = F_x;
        fx{i} = h_x;
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

    % Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + V*theta, theta in R^m
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
        F_u*K, F_u*(PSI - K*LAMBDA); ...
        zeros(length_Fu, n) F_u*PSI; ...
        F_x*(A+B*K) F_x*B*(PSI-K*LAMBDA)];
    h_w = [...
        h_x; ...
        h_x - F_x*LAMBDA_0; ...
        h_u - F_u*(PSI_0 - K*LAMBDA_0); ...
        h_u - F_u*PSI_0; ...
        h_x - F_x*B*(PSI_0-K*LAMBDA_0)];
    
    % Subtract out points due to disturbance (F_g,h_g)
%     [F_w_N0, h_w_N0] = pdiff(F_w, h_w, ...
%         [F_g zeros(length_Fg,m); ...
%         zeros(m, n) eye(m); ...
%         zeros(m, n) -eye(m)], ...
%         [h_g; ...
%         zeros(2*m,1)]);

    % Simplify the constraints
%     term_poly = polytope(F_w_N0, h_w_N0);
    term_poly = polytope(F_w, h_w);
    [F_w_N, h_w_N] = double(term_poly);
%     term_poly = Polyhedron(F_w_N0, h_w_N0); 
%     F_w_N = term_poly.A; % Inequality description { x | H*[x; -1] <= 0 }   
%     h_w_N = term_poly.b; % Inequality description { x | A*x <= b }

    disp('Terminal set Polyhedron:');
    term_poly

    % x-theta constraints:
    F_xTheta = F_w_N(:, 1:n);
    F_theta = F_w_N(:,n+1:n+m);
    f_xTheta = h_w_N;


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
    
        %%
    %==========================================================================
    % Compute uncertainty bounds
    %==========================================================================

    disp('Generating uncertainty bounds...');
    options = optimset('Display', 'final');
    uncertainty_upper_bounds = zeros(n, n + m + 1);
    uncertainty_lower_bounds = zeros(n, n + m + 1);
    
    X_vertices = extreme(polytope(blkdiag(F_x, F_u), [h_x; h_u]));
    length_X_vertices = size(X_vertices,1);
%     poly = Polyhedron(blkdiag(F_x, F_u), [h_x; h_u]);
%     X_vertices = poly.extreme();
%     length_X_vertices = size(X_vertices,1);
    X_vertices = [X_vertices ones(length_X_vertices,1)];

    for ind = 1:n
        disp(['sub_block ind: ', num2str(ind)]);

        if (sub_block(ind) > 0)
            s_vec = zeros(n,1); s_vec(ind) = 1;
            F_g_s = [1; -1]; h_g_s = [s_vec'*linprog(-s_vec, F_g, h_g); -s_vec'*linprog(s_vec, F_g, h_g, [], [], [], [], [], options)];
            F_delta_A = kron(F_g_s, X_vertices(:, logical(uncertainty_block(ind, :))));
            h_delta_A = repmat(h_g_s, length_X_vertices, 1);
            tempPol = Polyhedron(F_delta_A, h_delta_A);
            F_delta_A = tempPol.A; 
            h_delta_A = tempPol.b;

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
end
