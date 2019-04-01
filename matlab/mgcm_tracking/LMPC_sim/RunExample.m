%% TRACKING piecewise constant REFERENCE MPC example
close all;
clearvars;

%% Parameters
% Horizon length
N=50;
% Simulation length (iterations)
iterations = 1000;

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
% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;

umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];


%%
%  Shift the constraints for the linearised model for the value of the
%  working point
x_w = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
r0 = x_w(3);

F_u = [eye(m); -eye(m)]; h_u = [umax-r0; -umin+r0];
F_x = [eye(n); -eye(n)]; h_x = [xmax-x_w; -xmin+x_w];

% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);

run_F = [F_x zeros(length_Fx, m);...
        zeros(length_Fu,n) F_u];
run_h = [h_x;h_u];
%==========================================================================
% Compute maximally invariant set
%==========================================================================
L=PSI - K*LAMBDA;
disp('Computing and simplifying terminal set...');
F_w = [F_x zeros(length_Fx, m);
    zeros(length_Fx, n) F_x*LAMBDA; ...
    F_u*K, F_u*(L); ...
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
%%%%%%% NEW MPIS %%%%%%%
Ak=[A+B*K B*L; zeros(m,n) eye(m)];
Xc = Polyhedron(F_w,h_w);
term_poly2 = compute_MPIS(Xc,Ak);
MAI2=projection(term_poly2,1:2); % Maximal Admissible Invariant set projected on X
figure;
plot(MAI2);
F_w_N = term_poly2.A; % Inequality description { x | H*[x; -1] <= 0 }   
h_w_N = term_poly2.b; % Inequality description { x | A*x <= b }\
%%%%%%%%%%%%%%%%%%%%%%%
%%
% Kstable=[+3.0741 2.0957 0.1197 -0.0090]; % K stabilising gain from the papers
Kstable=-[3.0742   -2.0958   -0.1194    0.0089];

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];
%% Start simulation
sysHistory = [x_eq_init;u0(1:m,1)];
art_refHistory =  0;
true_refHistory = x_eq_ref;
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
x_init_true=x_eq_init+x_w; % init true sys state
x_ref_true=x_eq_ref+x_w;

x = x_w+x_eq_init; % real system input

tic;
for k = 1:(iterations)      
    %%
    fprintf('iteration no. %d/%d \n',k,iterations);
    % To give the values to the nominal mdel w.r.t. the point around whiuch
    % it is linearised around
    if k==1 
        x_eq = x_eq_init; 
    else
        x_eq = x - x_w;
    end
    COSTFUN = @(var) costFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x_eq,x_eq_ref,N,reshape(var(1:m),m,1),Q,R,P,T,K,LAMBDA,PSI);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x_eq,N,K,LAMBDA,PSI,F_x,h_x,F_u,h_u,F_w_N,h_w_N);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    
    % Implement first optimal control move and update plant states.
    [x, u] = getTransitionsTrue(x,c,x_w,r0,Kstable);
    
    % shift the output so that it's from the working point perspective
    % setpoint being [0;0;0;0]
    his = [x-x_w; u-r0]; % c: decision var; u-r0: delta u;
    
    % Save plant states for display.
    sysHistory = [sysHistory his]; 
    art_refHistory = [art_refHistory art_ref(1:m)];
    true_refHistory = [true_refHistory x_eq_ref];
    
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
legend({'input variable'},'Location','northeast')
grid on
xlabel('iterations');
ylabel('u');
title('Sys input');

% figure;
% plot_refs=plot(0:iterations, sysHistory(1:4,:), 'Linewidth',1.5);
% grid on;
% xlabel('iterations');
% ylabel('responses');
% plot_refs(1).Color = 'Yellow';
% plot_refs(2).Color = 'Blue';
% plot_refs(3).Color = 'Red';
% plot_refs(4).Color = 'Green';
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

