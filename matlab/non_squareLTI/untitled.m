%% Non-square LTI system 
A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
% for stady state parametrisation
Mtheta = [1, 0, 0, 0; 0, 1, 1, -2];
Ntheta =[1, 0];
n = size(A,1);
m = size(B,2);
Q = eye(n);
R = eye(m);
P = dare(A,B,Q,R);
T = 100*P;

%% Configuration
mgcm_bin_fname = 'nonsq.bin';

N_values = [3 ... % horizons for which to calculate the discriminating kernel (maximum disturbance invarinat set)
    % 60 120 240 ...
    ];

for N=N_values
    nonsq_bin_fname = ['nonsq_N' num2str(N) '.bin'];
    %% Get system definition:
    [A,B,C,d_0,Mtheta, Ntheta] = define_system();
    % Count number of states n, number of inputs m, number of outputs o:
    n = size(A,1);
    m = size(B,2);
    o = size(C,1);
    %% Get design parameters:
    [Q, R, ...
    dlqr_controlweight, maxadm_controlweight, ...
    mflow_max, mflow_min, ...
    prise_max, prise_min, ...
    throttle_max, throttle_min, ...
    throttle_rate_max, throttle_rate_min, ...
    u_max, u_min, ...
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
%     K = -dlqr(A, B, Q, dlqr_controlweight*R);

    % 'baseline' stabilizing feedback law
    K=[-3.0741 2.0957 0.1197 -0.0090]; %nominal feedback matrix from the LBMPC paper
%     Kt=K;
    % Terminal cost chosen as solution to DARE
    P = dare(A+B*K, B, Q, dlqr_controlweight*R);
    %%
    %==========================================================================
    % Define polytopic constraints on input F_u*x <= h_u and
    % state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
    %==========================================================================

    u_max;
    u_min;

    F_u = [eye(m); -eye(m)]; h_u = [u_max; -u_min];
    temp_max = [mflow_max; prise_max; throttle_max; throttle_rate_max];
    temp_min = [mflow_min; prise_min; throttle_min; throttle_rate_min];
%     temp_min =[-mflow_max; -prise_max; -throttle_max; -throttle_rate_max];
    F_x = [eye(n); -eye(n)]; h_x = [temp_max; -temp_min];
    F_g = [eye(n); -eye(n)]; h_g = [state_uncert; state_uncert]; % uncertainty Polyhedron

    % count the length of the constraints on input, states, and uncertainty:
    length_Fu = length(h_u);
    length_Fx = length(h_x);
    length_Fg = length(h_g);

    %% State constraints:
    % Points that after (A+B*K_t) get to (F_x,h_x) \ominus (F_g,h_g)
    [F_x_g, h_x_g] = double(polytope(F_x, h_x) - polytope(F_g, h_g));
%     tempPoly = Polyhedron(F_x, h_x) - Polyhedron(F_g, h_g);
%     F_x_g = tempPoly.A; % Inequality description { x | H*[x; -1] <= 0 }
%     h_x_g = tempPoly.b; % Inequality description { x | A*x <= b }
    
    Fx{1} = F_x;
    fx{1} = h_x;
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
