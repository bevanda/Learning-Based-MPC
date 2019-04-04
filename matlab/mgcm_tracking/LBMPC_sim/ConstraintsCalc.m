%% TRACKING piecewise constant REFERENCE MPC example
close all;
clearvars;

%% Parameters
% Horizon length
N=20;
% Simulation length (iterations)
iterations = 600;

%% Discrete time nominal model of the non-square LTI system for tracking
A = [1.01126321746508,-0.0100340214950357,6.46038913508018e-05,1.93716902346107e-07; ...
    0.0100340214950357,0.995515380253533,-0.0127681799951143,-5.57226765949308e-05; ...
    0,0,0.957038195891878,0.00792982548734094; ...
    0,0,-7.92982548734093,0.602405619103784];
B = [4.95338239742896e-07; ...
    -0.000193159646826652; ...
    0.0429618041081219; ...
    7.92982548734093];
C = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
% D = [0;0;0;0];
n = size(A,1);
m = size(B,2);
o = size(C,1);

% The initial conditions
x_eq_init = [-0.35;...
    -0.4;...
    0.0;...
    0.0];
%setpoint
x_eq_ref = [0.0;...
      0.0;...
      0.0;...
      0.0];

%%
%==========================================================================
% Generate steady-state parameterization
%==========================================================================
% MN = [Mtheta; 1, 0];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);
%%%%%%%%%%%%%%%%%%%%%
d_0 = [0,0,0,0]'; % inital disturbance guess
% Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + V*theta, theta in R^m
V_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = V_0(1:n);
PSI_0 = V_0(n+1:n+m);
%%%%%%%%%%%%%%%%%%%%%

%%
%==========================================================================
% Define a nominal feedback policy K and corresponding terminal cost
% 'baseline' stabilizing feedback law
Q = eye(n);
R = eye(m);

K = -dlqr(A, B, Q, R);
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
% terminal steady state cost
T = 1000; 


%%
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================
% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;

umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];
%%%%%%%%%%%%%%%%%%%%%%%
% To be calculated with the Taylor remainder theorem
% state_uncert = [0.0;0.0;0.0;0.0]; % just for testing
state_uncert = [0.02;5e-04;0;0]; % from max lin error from TaylorRemainder (Lagrange Error Bound)

%%%%%%%%%%%%%%%%%%%%%%%
%%
%  Shift the constraints for the linearised model for the value of the
%  working point
x_w = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
r0 = x_w(3);

% Shift the abs system constraints w.r.t. to the linearisation point
F_u = [eye(m); -eye(m)]; h_u = [umax-r0; -umin+r0];
F_x = [eye(n); -eye(n)]; h_x = [xmax-x_w; -xmin+x_w];
F_g = [eye(n); -eye(n)]; h_g = [state_uncert; state_uncert]; % uncertainty polytope
% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);
length_Fg = length(h_g);
% run_F = [F_x zeros(length_Fx, m);...
%         zeros(length_Fu,n) F_u];
% run_h = [h_x;h_u];
%% State constraints
temp = polytope(F_x, h_x) - polytope(F_g, h_g);
    [F_x_g, h_x_g] = double(temp);
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
% Compute maximal invariant set
%==========================================================================

% Terminal feedback policy for terminal set computations
maxadm_controlweight = 1/(umax^2); % r_i as the inverse of the square of the maximum permissible value for the corresponding u_i
K_t = -dlqr(A, B, Q, maxadm_controlweight*R);
%lambda=0.99; % λ ∈ (0, 1), λ can be chosen arbitrarily close to 1, the obtained
% invariant set can be used as a reliable polyhedral approximation to the maximal invariant set 
disp('Computing and simplifying terminal set...');
% extended state constraints

F_w = [ F_x zeros(length_Fx, m);
        zeros(length_Fx, n) F_x*LAMBDA; ...
        F_u*K_t, F_u*(PSI - K_t*LAMBDA); ...
        zeros(length_Fu, n) F_u*PSI; ...
        F_x_g*(A+B*K_t) F_x_g*B*(PSI-K_t*LAMBDA)];
h_w = [ h_x; ...
        h_x - F_x*LAMBDA_0; ...
        h_u - F_u*(PSI_0 - K_t*LAMBDA_0); ...
        h_u - F_u*PSI_0; ...
        h_x_g - F_x_g*B*(PSI_0-K_t*LAMBDA_0)];

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % from LMPC    
% F_w = [F_x zeros(length_Fx, m);
%     zeros(length_Fx, n) F_x*LAMBDA; ...
%     F_u*K, F_u*(PSI - K*LAMBDA); ...
%     zeros(length_Fu, n) F_u*PSI];
% 
% lambda=0.99; % λ ∈ (0, 1), λ can be chosen arbitrarily close to 1, the obtained
% % invariant set can be used as a reliable polyhedral approximation to the maximal invariant set 
% h_w = [...
%     h_x; ...
%     (h_x - F_x*LAMBDA_0)*lambda; ...
%     h_u - F_u*(PSI_0 - K*LAMBDA_0); ...
%     (h_u - F_u*PSI_0)]*lambda;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disturbance constraints of the extended state 
F_g_w = [F_g zeros(length_Fg,m); ...
        zeros(m, n) eye(m); ...
        zeros(m, n) -eye(m)];
h_g_w = [h_g; ...
        zeros(2*m,1)];
    

% calculating the robust positively invariant set    
[F_w_N0, h_w_N0] = pdiff(F_w, h_w, F_g_w, h_g_w);

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_w_N0 = F_w; h_w_N0 = h_w; % from LMPC
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simplify the constraints
term_poly = polytope(F_w_N0, h_w_N0);
[F_w_N, h_w_N] = double(term_poly);
% term_poly = Polyhedron(F_w_N01, h_w_N01);
% term_poly.minHRep(); % simplifying the polyhedron/constraints
% F_w_N1 = term_poly.A; h_w_N1 = term_poly.b;
disp('Terminal set Polyhedron:');
term_poly;
termProj = projection(term_poly,1:2);
% figure;
% plot(termProj);
disp('Terminal set Polyhedron 2:');
term_poly2 = polytope(F_w,h_w);
termProj2 = projection(term_poly2,1:2);
% figure;
% plot(termProj2);
disp('Terminal set Polyhedron diff:');
term_dif = term_poly-term_poly2;
isempty(term_dif);
termProjd = projection(term_dif,1:2);
% figure;
% plot(termProjd);
%% Calculation the n-reachable set
state_uncert = [0.02;5e-04;0;0]; % from max lin error from TaylorRemainder (Lagrange Error Bound)
F = [eye(n); -eye(n)]; h = [state_uncert; state_uncert];
W=Polyhedron(F,h);
W2=projection(W,1:2);
% figure;

% plot(W2);
AK=A+B*K;
% WR_10=fwd_reach_set(AK(1:2,1:2),W2,5);
% figure;
% plot(WR_10);
%% Full calc for RPI

% WR=reach_set_old(AK,W,5);
% figure;
% plot(WR.projection(1:2));
WR2=reach_set(AK,W,15);
figure;
plot(WR2.projection(1:2));

%%
sys = LTISystem('A', AK(1:2,1:2));

%% RPI calc
% WR_10=reach_set(AK(1:2,1:2),W2,10);
% figure;
% plot(WR_10);
% WR_20=reach_set_old(AK(1:2,1:2),W2,10);
% figure;
% plot(WR_20);

% 
% %% MPT fwd reach calc
% iset = sys.reachableSet('X', W2,'N',100,'direction','backward');
% figure;
% plot(iset);
% %% ROA calc
% WR_20=reach(AK(1:2,1:2),W2,10);
% figure;
% plot(WR_20);
% 
% %% MPT back reach calc
% WR_20=reach(AK(1:2,1:2),W2,10);
% figure;
% plot(WR_20);
%%
% %==========================================================================
% % Generate inequality constraints
% %==========================================================================
% 
% length_Fw = size(F_w_N, 1);
% 
% Aineq = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, N*m+m);
% bineq = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, 1);
% b_crx = zeros((N-1)*length_Fx+N*length_Fu+length_Fw, n);
% % help variables
% L_i = zeros(n, N*m); % width of the state tube R_i 
% KL_i = zeros(m, N*m); % width of the input tube KR_i
% disp('Generating constraints on inputs...');
% d_i = zeros(n,1);
% for ind = 1:N
%     disp(['u ind: ', num2str(ind)]);
% 
%     KL_i = K*L_i;
%     KL_i(:, (ind-1)*m + (1:m)) = eye(m);
% 
%     Aineq((ind-1)*length_Fu + (1:length_Fu),1:N*m) = F_u*KL_i;
%     bineq((ind-1)*length_Fu + (1:length_Fu)) = h_u - F_u*K*d_i;
%     b_crx((ind-1)*length_Fu + (1:length_Fu),:) = -F_u*K*(A+B*K)^(ind-1);
% 
%     L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];
%     d_i = (A+B*K)*d_i + d_0;
% end
% 
% L_i = zeros(n, N*m);
% disp('Generating constraints on states...');
% d_i = d_0;
% for ind = 1:N
%     disp(['x ind: ', num2str(ind)]);
%     L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];
% 
%     if ind == 1
%         disp('Generating terminal constraints on states...');
%         Aineq(N*length_Fu + (1:length_Fw), :) = F_w_N*[L_i zeros(n,m); zeros(m,N*m) eye(m)];
%         bineq(N*length_Fu + (1:length_Fw)) = h_w_N - F_w_N*[d_i; zeros(m,1)];
%         b_crx(N*length_Fu + (1:length_Fw),:) = -F_w_N*[(A+B*K)^(ind); zeros(m,n)];
% 
%     else
% 
%         Aineq(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx),1:N*m) = F_x*L_i;
%         bineq(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx)) = h_x - F_x*d_i;
%         b_crx(length_Fw + N*length_Fu + (ind-2)*length_Fx + (1:length_Fx),:) = -F_x*(A+B*K)^(ind);
% 
%     end
% 
%     d_i = (A+B*K)*d_i + d_0;
% end
% 
% ind = N;
% L_i = [(A+B*K)^(ind-1)*B L_i(:, 1:(N-1)*m)];
% 
% CONSTRAINT_COUNT = length(bineq);
% disp('Removing redundant inequality constraints...');
% temp_tope = polytope(Aineq, bineq);
% [Aineq, bineq] = double(temp_tope);