%% TRACKING piecewise constant REFERENCE MPC example
close all;
clearvars;

%% Parameters
% Set the prediction horizon:
N = 50;
% Simulation length (iterations)
iterations = 600;
% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;
%% Closed-Loop Simulation
% working point 
xw = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
% The initial conditions
x = [-0.35;...
    -0.4;...
    0.0;...
    0.0];
%setpoint
xs = [0.0;...
      0.0;...
      0.0
      0.0];

options = optimoptions('fmincon','Algorithm','sqp','Display','final');

% LB =[];
% UB = [];

A = [1.01125000000000,0.0100000000000000,0,0;...
    0.0100000000000000,0.995555557627778,-0.0129903810567666,0;...
    0,0,1,0.0100000000000000;...
    0,0,-10,0.552786404500042];
B = [0;0;0;10];
C = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
D = [0;0;0;0];

n = size(A,1);
m = size(B,2);
o = size(C,1);


% MN = [Mtheta; 1, 0];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);
%%%%%%%%%%%%%%%%%%%%%
d_0 = [0,0,0,0]';
% Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + V*theta, theta in R^m
V_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = V_0(1:n);
PSI_0 = V_0(n+1:n+m);
%%%%%%%%%%%%%%%%%%%%%
Q = eye(n);
R = eye(m);
%%
%==========================================================================
% Define a nominal feedback policy K and corresponding terminal cost
% 'baseline' stabilizing feedback law
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

umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];


F_u = [eye(m); -eye(m)]; h_u = [umax; -umin];
% deduce the working point from the constraints for the linearised model
F_x = [eye(n); -eye(n)]; h_x = [xmax-xw; -(xmin-xw)];

% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);

run_F = [F_x zeros(length_Fx, m);...
        zeros(length_Fu,n) F_u];
run_h = [h_x;h_u];
%%
%==========================================================================
% Compute maximally invariant set
%==========================================================================

disp('Computing and simplifying terminal set...');
F_w = [F_x zeros(length_Fx, m);
    zeros(length_Fx, n) F_x*LAMBDA; ...
    F_u*K, F_u*(PSI - K*LAMBDA); ...
    zeros(length_Fu, n) F_u*PSI];

lambda=0.99; % λ ∈ (0, 1), λ can be chosen arbitrarily close to 1, the obtained
% invariant set can be used as a reliable polyhedral approximation to the maximal invariant set 
h_w = [...
    h_x; ...
    (h_x - F_x*LAMBDA_0)*lambda; ...
    h_u - F_u*(PSI_0 - K*LAMBDA_0); ...
    (h_u - F_u*PSI_0)]*lambda;

F_w_N0 = F_w; h_w_N0 = h_w;

% Simplify the constraints
term_poly = polytope(F_w_N0, h_w_N0);
[F_w_N, h_w_N] = double(term_poly);
%     term_poly = Polyhedron(F_w_N0, h_w_N0); 
%     F_w_N = term_poly.A; % Inequality description { x | H*[x; -1] <= 0 }   
%     h_w_N = term_poly.b; % Inequality description { x | A*x <= b }
disp('Terminal set Polyhedron:');
term_poly
% MAI=projection(term_poly,1:2); % Maximal Admissible Invariant set rpojected on X
% plot(MAI);
% x-theta constraints:
% F_xTheta = F_w_N;
% F_x = F_w_N(:, 1:n);
% F_theta = F_w_N(:,n+1:n+m);
% f_xTheta = h_w_N;
LB = [ones(m*N,1)*u_min; ones(m,1)*(-inf)];
UB = [ones(m*N,1)*u_max; ones(m,1)*(+inf)];
u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];
%% Start simulation
sysHistory = [x;u0(1:m,1)];
art_refHistory =  0;
true_refHistory = xs;
tic;
for k = 1:(iterations)
    COSTFUN = @(var) costFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,xs,N,reshape(var(1:m),m,1),Q,R,P,T,K,LAMBDA,PSI);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,N,K,LAMBDA,PSI,run_F,run_h,F_w_N,h_w_N);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    % Implement first optimal control move and update plant states.
    x= getTransitions(x, c, K); %-K*(x-xs)+u_opt
    his=[x; c]; % c = delta u
    % Save plant states for display.
    sysHistory = [sysHistory his]; 
    art_refHistory = [art_refHistory art_ref(1:m)];
    true_refHistory = [true_refHistory xs];
end
toc
%% Plot

figure;
subplot(n+m,1,1);
plot(0:iterations,sysHistory(1,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x1');
title('mass flow');
subplot(n+m,1,2);
plot(0:iterations,sysHistory(2,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x2');
title('pressure rise');
subplot(n+m,1,3);
plot(0:iterations,sysHistory(3,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x3');
title('throttle');
subplot(n+m,1,4);
plot(0:iterations,sysHistory(4,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x4');
title('throttle rate');
subplot(n+m,1,5);
plot(0:iterations,sysHistory(5,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('u');
title('delta u input from MPC');

figure;
plot_refs=plot(0:iterations,art_refHistory(1,:), 0:iterations, true_refHistory(1,:),0:iterations,sysHistory(1,:),'Linewidth',1.5);
grid on
xlabel('iterations');
% ylabel('references');
title('Artificial vs true reference vs state response');
legend({'artifical reference','real reference', 'mass flow response'},'Location','northeast')
plot_refs(1).LineStyle='--';
plot_refs(2).LineStyle='-.';
plot_refs(1).Color='Red';
plot_refs(2).Color='Black';
plot_refs(3).Color='Blue';

figure;
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1,'Marker','.');
grid on
xlabel('x1');
ylabel('x2');
title('State space');

