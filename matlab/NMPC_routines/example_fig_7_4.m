% -------------------------------------------------------------------------
%                        Example: inverted pendulum
% -------------------------------------------------------------------------
function main
    clearvars;
% -------------------------------------------------------------------------
% Input parameter
% -------------------------------------------------------------------------

    dimension_u = 1;                                        % Dimension of u(t)
    dimension_x = 4;                                        % Dimension of x(t)

    g = 9.81;                                               % Gravitation
    k = 0.1;                                                % Friction

    A_continuous = [0 1 0 0; g -k 0 0; 0 0 0 1; 0 0 0 0];   % Coefficient matrix A
    B_continuous = [0; 1; 0; 1];                            % Coefficient matrix B

    Q = 2.0*eye(dimension_x, dimension_x);                  % Weighting matrix Q
    R = 4.0*eye(dimension_u, dimension_u);                  % Weighting matrix R
    endweight = 1;                                          % Weighting of endpoint

    N_u = 9;                                                % Length of the control horizon
    N_x = N_u;                                              % Length of the prediction horizon

    t_akt = 0;                                              % Actual time
    t_s = 0;                                                % Time to start MPC

    x1 = 0.0;                                               % Initial values
    x2 = 0.0;
    x3 = -4.0;
    x4 = -1.0;

    grid_diameter = 0.2;                                    % Diameter of grid for initial values
                                                            %   in norm infty
    gridnodes_per_dimension = 5;                            % Number of gridnodes per dimension
                                                            %   equally distributed

    T = 0.5;                                                % Sampling period

    iplot = 0;                                              % Switch for plotting
                                                            %  1-plotting
                                                            %  0-no plots
    N_plot = 19;                                            % Application horizon
    Nmvals = [1 2 3 4 5 6 7 8 9];                           % Values for m-step feedback

% -------------------------------------------------------------------------
% Setting up grid
% -------------------------------------------------------------------------
    x_t_akt = [];
    [x1_grid, x2_grid, x3_grid, x4_grid] = ...
        ndgrid(-grid_diameter/2:grid_diameter/(gridnodes_per_dimension-1):grid_diameter/2, ...
                -grid_diameter/2:grid_diameter/(gridnodes_per_dimension-1):grid_diameter/2, ...
                -grid_diameter/2:grid_diameter/(gridnodes_per_dimension-1):grid_diameter/2, ...
                -grid_diameter/2:grid_diameter/(gridnodes_per_dimension-1):grid_diameter/2);
    for i = 1:gridnodes_per_dimension
        for j = 1:gridnodes_per_dimension
            for k = 1:gridnodes_per_dimension
                for l = 1:gridnodes_per_dimension
                    x_t_akt = [[x_t_akt]; [x1 x2 x3 x4]+[x1_grid(i,j,k,l) x2_grid(i,j,k,l) x3_grid(i,j,k,l) x4_grid(i,j,k,l)]];
                end
            end
        end
    end
    x_t_akt = x_t_akt';

    for i = 1:size(x_t_akt,2)
        fprintf('Computing trajectory %d of %d\n',i,size(x_t_akt,2));
% -------------------------------------------------------------------------
% Solving linear MPC problem
% -------------------------------------------------------------------------
        for N_m=Nmvals
            [alpha_sequence, V_sequence, l_sequence] = ...
                    programm_AB(A_continuous, B_continuous, T, Q, R, endweight, ...
                                N_u, N_x, N_m, N_plot, t_akt, t_s, x_t_akt(:,i),...
                                dimension_u, dimension_x, iplot);

            alpha_sum(i,N_m) = (V_sequence(1)-V_sequence(length(alpha_sequence)+1)) ...
                           /sum(l_sequence);
            alpha_min(i,N_m) = min(alpha_sequence);
            alpha_mean(i,N_m) = sum(alpha_sequence)/length(alpha_sequence);

            V_sequence_1(N_m) = V_sequence(1);
            V_sequence_last(N_m) = V_sequence(length(alpha_sequence)+1);
            l_sum(N_m) = sum(l_sequence);

            alpha_number(N_m) = length(alpha_sequence);
        end
% -------------------------------------------------------------------------
% Output of characteristic values for suboptimality
% -------------------------------------------------------------------------
        for N_m=Nmvals
            for j = 1:dimension_x
                fprintf('x_%i: %f  ', ...
                         j, x_t_akt(j,i));
            end
            fprintf('N_m=%d:  %f  (%d)\n',...
                     N_m, alpha_sum(i,N_m), N_m*alpha_number(N_m));
%                     N_m, alpha_min(i,N_m), N_m*alpha_number(N_m));
%                     N_m, alpha_mean(i,N_m), N_m*alpha_number(N_m));
        end
        if (i == 1)
            figure(1);
            clf;
            hold on;
%            plot(Nmvals, alpha_sum(i,Nmvals), '--x','Linewidth',1)
            plot(Nmvals, alpha_min(i,Nmvals), '--x','Linewidth',1)
%            plot(Nmvals, alpha_mean(i,Nmvals), '--x','Linewidth',1)
            title('\fontsize{14} Suboptimality estimate along the trajectory');
            xlabel('control horizon m');
            ylabel('suboptimality degree \alpha');
            axis([1 max(Nmvals) 0 1]);
        else
            figure(1);
%            plot(Nmvals, alpha_sum(i,Nmvals), '--x','Linewidth',1)
            plot(Nmvals, alpha_min(i,Nmvals), '--x','Linewidth',1)
%            plot(Nmvals, alpha_mean(i,Nmvals), '--x','Linewidth',1)
        end
    end
    figure(1);
%    plot(Nmvals, min(alpha_sum(:,Nmvals)), 'g-x','Linewidth',2,'MarkerFaceColor','r')
    plot(Nmvals, min(alpha_min(:,Nmvals)), 'g-x','Linewidth',2,'MarkerFaceColor','r')
%    plot(Nmvals, min(alpha_mean(:,Nmvals)), 'g-x','Linewidth',2,'MarkerFaceColor','r')
%    plot(Nmvals, sum(alpha_sum(:,Nmvals))/size(alpha_sum,1), 'r-x','Linewidth',2,'MarkerFaceColor','r')
%    plot(Nmvals, sum(alpha_min(:,Nmvals))/size(alpha_min,1), 'r-x','Linewidth',2,'MarkerFaceColor','r')
%    plot(Nmvals, sum(alpha_mean(:,Nmvals))/size(alpha_mean,1), 'r-x','Linewidth',2,'MarkerFaceColor','r')
    figure(2);
    clf;
%    plot(Nmvals, min(alpha_sum(:,Nmvals)), 'g-x','Linewidth',2,'MarkerFaceColor','r')
    plot(Nmvals, min(alpha_min(:,Nmvals)), 'g-x','Linewidth',2,'MarkerFaceColor','r')
%    plot(Nmvals, min(alpha_mean(:,Nmvals)), 'g-x','Linewidth',2,'MarkerFaceColor','r')
    title('\fontsize{14} Suboptimality estimate for solutions eminating a given set');
    xlabel('control horizon m');
    ylabel('suboptimality degree \alpha');
    axis([1 max(Nmvals) 0 1]);



% -------------------------------------------------------------------------
% Plot of suboptimality estimates for different m-step feedbacks
% -------------------------------------------------------------------------
%     figure(1);
%         clf;
%         hold on;
%         plot(Nmvals, alpha_gesamt(Nmvals), 'b--o','Linewidth',2,'MarkerFaceColor','b');
%         plot(Nmvals, alpha_mean(Nmvals), 'r--s','Linewidth',2,'MarkerFaceColor','r');
%         plot(Nmvals, alpha_min(Nmvals), 'g--x','Linewidth',2,'MarkerFaceColor','g');
%         title('\fontsize{12} Suboptimality estimate along the trajectory');
%         legend('Overall \alpha', 'Mean of \alpha', 'Minimum of \alpha');
%         xlabel('m');
%         ylabel('\alpha');
%         axis([1 10 0.5 1]);
%         hold off;
%
%     figure(2);
%         clf;
%         plot(Nmvals,alpha(Nmvals),'k--o');
%         title('\fontsize{12} Suboptimality estimate for the first step');
%         xlabel('m');
%         ylabel('\alpha');
%         hold off;








% -------------------------------------------------------------------------
%           Model predictive control using linear programming
% -------------------------------------------------------------------------

% Requirement: linear continuous-time control system

function [alpha_sequence,V_sequence,l_sequence] = ...
            programm_AB(A_continuous, B_continuous, T, Q, R, endweight, ...
                        N_u, N_x, N_m, N_plot, t_akt, t_s, x_t_akt, ...
                        dimension_u, dimension_x, iplot)

% -------------------------------------------------------------------------
% Discretization of the problem
% -------------------------------------------------------------------------
    A = c2d_A(A_continuous, B_continuous, T);
    B = c2d_B(A_continuous, B_continuous, T);

% -------------------------------------------------------------------------
% Variables for graphical output
% -------------------------------------------------------------------------
    T_sampling = [t_s : t_s+N_plot]*T;
    u_MPC = [zeros(dimension_u, N_plot-1)];
    x_MPC = [zeros(dimension_x, N_plot)];
    t_continuous  = [ ];
    x_continuous  = [ ];

% -------------------------------------------------------------------------
%  Setup of the linear programm to solve the MPC problem
% -------------------------------------------------------------------------
%    options = optimset('LargeScale', 'off', 'Simplex', 'on');
%     options = optimset('LargeScale', 'off', 'dual-simplex', 'off','display','off');
   options = optimset('LargeScale', 'on');

    f = vector_f(N_u,N_x,dimension_u,dimension_x,endweight);

    E1 = matrix_E1(A,B,Q,N_u,N_x,dimension_u,dimension_x);
    E2 = matrix_E2(R,N_u,dimension_u);
    E3 = matrix_E3(N_x,dimension_x);
    E4 = matrix_E4(N_u,dimension_u);
    Z1 = matrix_Z1(N_u,N_x,dimension_u,dimension_x);
    Z2 = matrix_Z2(N_u,N_x,dimension_u,dimension_x);

    E = [E1 E3 Z2; -E1 E3 Z2; E2 Z1 E4; -E2 Z1 E4];

% -------------------------------------------------------------------------
%  Variables for suboptimality estimation
% -------------------------------------------------------------------------
    alpha_sequence =[ ];
    V_sequence = [ ];
    V_previous = 0;
    V_actual = 0;
    l_sequence = [ ];
    l_previous = 0;
    l_actual = 0;

% -------------------------------------------------------------------------
%  MPC routine
% -------------------------------------------------------------------------
    t = t_s;
    sampling_step = 1;
    x_MPC(:, sampling_step) = x_t_akt;
    while t <= N_plot-1
        V_previous = V_actual;
        l_previous = l_actual;

        d1 = vector_d1(A, Q, N_x, dimension_x, t_s, t_akt, x_t_akt);
        d2 = vector_d2(N_u, dimension_u);
        d = [d1; -d1; d2; d2];
% -------------------------------------------------------------------------
%   Solving the linear programm
% -------------------------------------------------------------------------
        z = linprog(f, E, d, [ ], [ ], [ ], [ ], [ ], options);
        u = [z(1:dimension_u*N_u,:); zeros(dimension_u,1)];
%        min = f*z;
%        V_actual_opt = min+norm(Q*x_t_akt,1);
% -------------------------------------------------------------------------
%   Calculating optimal value function value on the control interval and
%   stage cost value on the implementation interval
% -------------------------------------------------------------------------
        if (iplot ==1)
            V_actual = evalV(A, B, Q, R, N_u+1, dimension_u, x_t_akt, u, ...
                             endweight)
        else
            V_actual = evalV(A, B, Q, R, N_u+1, dimension_u, x_t_akt, u, ...
                             endweight);
        end
        l_actual = evalV(A, B, Q, R, N_m, dimension_u, x_t_akt, u, 1);
% -------------------------------------------------------------------------
%   Implementation of m-step MPC feedback
% -------------------------------------------------------------------------
        for step_counter = 1:min([N_m; (N_plot-t)])
% -------------------------------------------------------------------------
%       Simulating the discrete time control system using MPC feedback
% -------------------------------------------------------------------------
            x  = A*x_t_akt+B*u((step_counter-1)*dimension_u+...
                    1:step_counter*dimension_u,:);
            x_MPC(:,sampling_step+1) = x;
            u_MPC(:,sampling_step) = u((step_counter-1)*dimension_u+...
                    1:step_counter*dimension_u,:);
% -------------------------------------------------------------------------
%       Simulating the continuous time control system using MPC feedback
% -------------------------------------------------------------------------
            [t_continuousl,x_continuousl] = simulate_cont(A_continuous, ...
                    B_continuous, u((step_counter-1)*dimension_u+...
                    1:step_counter*dimension_u, :), x_t_akt, ...
                    [t*T, (t+1)*T]);
            t_continuous = [t_continuous; t_continuousl];
            x_continuous = [x_continuous; x_continuousl];
% -------------------------------------------------------------------------
%       Shift of the MPC problem
% -------------------------------------------------------------------------
            x_t_akt = x;
            t = t+1;
            sampling_step = sampling_step+1;
        end
% -------------------------------------------------------------------------
%   Computing alpha values
% -------------------------------------------------------------------------
        if (V_previous > 0)
            alpha_sequence = [alpha_sequence; ...
                              (V_previous-V_actual)/l_previous];
            V_sequence = [V_sequence; V_actual];
            l_sequence = [l_sequence; l_previous];
        else
            V_sequence = [V_actual];
        end
    end

% -------------------------------------------------------------------------
% Graphical output
% -------------------------------------------------------------------------
    if (iplot ==1)
        figure;
            plot(t_continuous, x_continuous);
            hold on;
            plot(T_sampling, x_MPC, '+');
            stairs(T_sampling(1:length(T_sampling)-1), u_MPC, 'b--');
    end

% -------------------------------------------------------------------------
% Evaluating cost functional
% -------------------------------------------------------------------------
function V = evalV(A, B, Q, R, N, dimension_u, x, u, endweight)
    V = 0;
    for k = 0:N-2
        V = V+norm(Q*x, 1)+...
            norm(R*u(k*dimension_u+1:(k+1)*dimension_u, :), 1);
        x  = A*x+B*u(k*dimension_u+1:(k+1)*dimension_u, :);
    end
    V = V+endweight*(norm(Q*x, 1)+...
        norm(R*u((N-1)*dimension_u+1:N*dimension_u, :), 1));

% -------------------------------------------------------------------------
% Compute matrix A of the discrete time control system
% -------------------------------------------------------------------------
function A = c2d_A(A_continuous,B_continuous,t_abt)
    [ma,na] = size(A_continuous);
    [mb,nb] = size(B_continuous);
    s = expm([[A_continuous B_continuous]*t_abt; zeros(nb,na+nb)]);
    A = s(1:na,1:na);

% -------------------------------------------------------------------------
% Compute matrix B of the discrete time control system
% -------------------------------------------------------------------------
function B = c2d_B(A_continuous,B_continuous,t_abt)
    [ma,na] = size(A_continuous);
    [mb,nb] = size(B_continuous);
    s = expm([[A_continuous B_continuous]*t_abt; zeros(nb,na+nb)]);
    B = s(1:na,na+1:na+nb);

% -------------------------------------------------------------------------
% Construction matrix E1
% -------------------------------------------------------------------------
function E1 = matrix_E1(A,B,Q,N_u,N_x,m,n)
    E1 = [zeros(n*N_u+n*(N_x-N_u),m*N_u)];
        E11 = matrix_E11(B,Q,m,n);
        E12_to_E1Nx = matrix_E12_to_E1Nx(A,B,Q,N_u,N_x,m,n);
        E1_1 = matrix_E1_1(N_u,N_x,E11,m,n);
        E1_2 = matrix_E1_2(N_u,N_x,E12_to_E1Nx,m,n);
    E1 = E1_1+E1_2;

% -------------------------------------------------------------------------
% Construction matrix E11
% -------------------------------------------------------------------------
function E11 = matrix_E11(B,Q,m,n)
    E11 = [zeros(n,m)];
    for i = 1:m
        for j = 1:n
            E11(j,i) = Q(j,j)*B(j,i);
        end
    end

% -------------------------------------------------------------------------
% Construction matrix E12_to_E1Nx
% -------------------------------------------------------------------------
function E12_to_E1Nx = matrix_E12_to_E1Nx(A,B,Q,N_u,N_x,m,n)
    E12_to_E1Nx = [zeros(n*(N_x-1),m)];
    s = zeros(1,n);
    for i = 1:m
        for j = 1:N_x-1
            A1 = A^j;
            for k = 1:n
                for l = 1:n
                    s(l) = A1(k,l)*B(l,i);
                end
                E12_to_E1Nx((j-1)*n+k,i) = Q(k,k)*sum(s);
            end
        end
    end

% -------------------------------------------------------------------------
% Construction matrix E1_1
% -------------------------------------------------------------------------
function E1_1 = matrix_E1_1(N_u,N_x,E11,m,n)
    E1_1 = [zeros(n*(N_u)+n*(N_x-N_u),m*N_u)];
    for i = 1:N_u
        for j = 1:m
            for k = 1:n
                E1_1((i-1)*n+k,(i-1)*m+j) = E11(k,j);
            end
        end
    end


% -------------------------------------------------------------------------
% Construction matrix E1_2
% -------------------------------------------------------------------------
function E1_2 = matrix_E1_2(N_u,N_x,E12_to_E1Nx,m,n)
    E1_2 = [zeros(n*N_x,m*N_u)];
    for i = 1:N_u
        for j = 1:m
            for k = 1:n*(N_x-i)
                E1_2(n*i+k,(i-1)*m+j) = E12_to_E1Nx(k,j);
            end
        end
    end

% -------------------------------------------------------------------------
% Construction matrix E2
% -------------------------------------------------------------------------
function E2 = matrix_E2(R,N_u,m)
    E2 = [zeros(m*N_u,m*N_u)];
    E21 = matrix_E21(R,m);
    for i = 1:N_u
        for j = 1:m
            for k = 1:m
                E2((i-1)*m+k,(i-1)*m+j) = E21(k,j);
            end
        end
    end

% -------------------------------------------------------------------------
% Construction matrix E21
% -------------------------------------------------------------------------
function E21 = matrix_E21(R,m)
    E21 = [zeros(m,m)];
    E21 = R;

% -------------------------------------------------------------------------
% Construction matrix E3
% -------------------------------------------------------------------------
function E3 = matrix_E3(N_x,n)
    E3 = [zeros(n*N_x,n*N_x)];
    E3 = diag(-ones(1,n*N_x));

% -------------------------------------------------------------------------
% Construction matrix E4
% -------------------------------------------------------------------------
function E4 = matrix_E4(N_u,m)
    E4 = [zeros(m*N_u,m*N_u)];
    E4 = diag(-ones(1,m*N_u));

% -------------------------------------------------------------------------
% Construction matrix Z1
% -------------------------------------------------------------------------
function Z1 = matrix_Z1(N_u,N_x,m,n)
    Z1 = [zeros(m*N_u,n*N_x)];

% -------------------------------------------------------------------------
% Constructing matrix Z2
% -------------------------------------------------------------------------
function Z2 = matrix_Z2(N_u,N_x,m,n)
    Z2 = [zeros(n*N_x,m*N_u)];

% -------------------------------------------------------------------------
% Constructing vector d1
% -------------------------------------------------------------------------
function d1 = vector_d1(A,Q,N_x,dimension_x,t_s,t_akt,x_t_akt)
    d1 = [zeros(dimension_x*N_x, 1)];
    for i = 1:N_x
        for j = 1:dimension_x
            A1 = A^i*A^(t_s-t_akt)*x_t_akt;
            d1((i-1)*dimension_x+j) =-Q(j, j)*A1(j,1);
        end
    end

% -------------------------------------------------------------------------
% Constructing vector d2
% -------------------------------------------------------------------------
function d2 = vector_d2(N_u,m)
    d2 = [zeros(m*N_u,1)];

% -------------------------------------------------------------------------
% Constructing vector f
% -------------------------------------------------------------------------
function f = vector_f(N_u,N_x,m,n,endweight)
    f = [zeros(1,2*m*N_u+n*N_x)];
        f1 = [zeros(1,m*N_u)];
        f2 = [ones(1,n*(N_x-1)) endweight*ones(1,n)];
        f3 = [ones(1,m*N_u)];
    f = [f1 f2 f3];

% -------------------------------------------------------------------------
% Simulating the continuous time system
% -------------------------------------------------------------------------
function [t, y] = simulate_cont(A_continuous, B_continuous, u, x0, tt)
    [t, y] = ode45(@(t, x)rhs(t, x, A_continuous, B_continuous, u), tt, x0);
function y = rhs(t, x, A_continuous, B_continuous, u)
    y = A_continuous*x+B_continuous*u;



