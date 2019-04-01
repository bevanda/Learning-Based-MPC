%% TRACKING piecewise constant REFERENCE MPC example
clearvars;
close all;

%% Parameters
% Set the prediction horizon:
N = 3;

% The initial conditions
x = [0;...
    -2];
%setpoint
xs = [4.95;...
      0.0];
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
% Simulation length (iterations)
iterations = 90;
% Optimization variable bounds not constraints per say
u_min =-0.3; u_max= 0.3;
LB = [ones(2*N,1)*u_min; ones(2,1)*(-inf)];
UB = [ones(2*N,1)*u_max; ones(2,1)*(+inf)];
% UB = [];
% LB = [];

A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
n = size(A,1);
m = size(B,2);
o = size(C,1);

Q = diag([1,1]);
R = diag([1,1]);

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];

% MN = [Mtheta; 1, 0];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);

%%%%%%%%%%%%%%%%%%%%%
% UNDER DISTURBANCE %
d_0 = [0,0]';
% Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + V*theta, theta in R^m
V_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = V_0(1:n);
PSI_0 = V_0(n+1:n+m);
%%%%%%%%%%%%%%%%%%%%%

%%
%==========================================================================
% Define a nominal feedback policy K and corresponding terminal cost

% 'baseline' stabilizing feedback law
K = -dlqr(A, B, Q, R);
u_min=-0.3;
max_admissible_ctrl_weight=1/(u_min^2);
K_t = -dlqr(A, B, Q, max_admissible_ctrl_weight*R);
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
% terminal steady state cost
T = 100*P;
%%
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================

umax = [0.3;0.3]; umin = [-0.3;-0.3];
xmax = [5; 5]; xmin = [-5; -5];

F_u = [eye(m); -eye(m)]; h_u = [umax; -umin];
F_x = [eye(n); -eye(n)]; h_x = [xmax; -xmin];

% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);

run_F = [F_x zeros(length_Fx,m);...
        zeros(length_Fu,n) F_u];
run_h = [h_x;h_u];
%%
%==========================================================================
% Compute maximally invariant set
%==========================================================================

disp('Computing and simplifying terminal set...');
L = (PSI - K*LAMBDA);
L0 = (PSI_0 - K*LAMBDA_0); % when being under inital disturbance
F_w = [F_x zeros(length_Fx, m);
    zeros(length_Fx, n) F_x*LAMBDA; ...
    F_u*K, F_u*L; ...
    zeros(length_Fu, n) F_u*PSI];
% contract the h values for the artificial steady state by a scalar λ ∈ (0, 1)
lambda=0.99;
% h_x = lambda*h_x;
% h_u = lambda*h_u;
h_w = [...
    h_x; ...
    lambda*(h_x - F_x*LAMBDA_0); ...
    h_u-F_u*L0; ...
    lambda*(h_u- F_u*PSI_0)];
F_w_N0 = F_w; h_w_N0 = h_w;

% Simplify the constraints
% term_poly = polytope(F_w_N0, h_w_N0);
% [F_w_N, h_w_N] = double(term_poly);
term_poly = Polyhedron(F_w_N0, h_w_N0); 
F_w_N = term_poly.A; % Inequality description { x | H*[x; -1] <= 0 }   
h_w_N = term_poly.b; % Inequality description { x | A*x <= b }\

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVEX POLYHAEDRON Wλ = {w = (x, θ) : (x, Kx + Lθ) ∈ Z, Mθθ ∈ λZ}.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Terminal set Polyhedron:');
term_poly
MAI=projection(term_poly,1:n); % Maximal Admissible Invariant set projected on X
plot(MAI);
% x-theta constraints:
F_xTheta = F_w_N;
F_x = F_w_N(:, 1:n);
F_theta = F_w_N(:,n+1:n+m);
f_xTheta = h_w_N;

%% Cost Calculation
% Start simulation
sysHistory = [x;u0(1:2,1)];
art_refHistory = LAMBDA*theta0;
true_refHistory = xs;

for k = 1:(iterations)
    xs = set_ref(k);
    
    COSTFUN = @(var) costFunction(reshape(var(1:end-2),2,N),reshape(var(end-1:end),2,1),x,xs,N,reshape(var(1:2),2,1),P,T,K,LAMBDA,PSI);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-2),2,N),reshape(var(end-1:end),2,1),x,N,K,LAMBDA,PSI,run_F,run_h,F_xTheta,f_xTheta);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-1:end),2,1);
    u = reshape(opt_var(1:2),2,1);
    art_ref = Mtheta*theta_opt;
    % Implement first optimal control move and update plant states.
    x= getTransitions(x, u); %-K*(x-xs)+u_opt
    his=[x; u];
    % Save plant states for display.
    sysHistory = [sysHistory his]; 
    art_refHistory = [art_refHistory art_ref(1:2)];
    true_refHistory = [true_refHistory xs];
end


%% Plot

figure;
subplot(3,1,1);
plot(0:iterations,sysHistory(1,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x1');
title('x1');
subplot(3,1,2);
plot(0:iterations,sysHistory(2,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x2');
title('x2');
subplot(3,1,3);
plot(0:iterations,sysHistory(3,:),0:iterations,sysHistory(4,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('u');
title('inputs');

figure;
plot_refs=plot(0:iterations,art_refHistory(1,:), 0:iterations, true_refHistory(1,:),0:iterations,sysHistory(1,:),'Linewidth',1.5);
grid on
xlabel('iterations');
% ylabel('references');
title('Artificial vs true reference vs state response');
legend({'artifical reference','real reference', 'state response'},'Location','northeast')
plot_refs(1).LineStyle='--';
plot_refs(2).LineStyle='-.';
plot_refs(1).Color='Red';
plot_refs(2).Color='Black';
plot_refs(3).Color='Blue';

figure;
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1.5,'Marker','.');
grid on
xlabel('x1');
ylabel('x2');
title('State space');

%% Helper functions

%set reference depending on the iteration
function [xs] = set_ref(ct)
    if ct <=30 
        xs=[4.95;0];
    elseif ct > 30 && ct <= 60
        xs=[-5.5;0];
    else
        xs=[2;0];
    end
end