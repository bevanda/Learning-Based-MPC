%% TRACKING piecewise constant REFERENCE MPC example
close all;
clearvars;

%% Parameters
N=10;
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
%  Shift the constraints for the linearised model for the value of the
%  working point
x_w = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
r0 = x_w(3);
% Kstable=[+3.0741 2.0957 0.1197 -0.0090]; % K stabilising gain from the papers
Kstable=[3.0742   -2.0958   -0.1194    0.0089];

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];
%% Start simulation
sysHistory = [x_eq_init;u0(1:m,1)];
true_refHistory = x_eq_ref;
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
x_init_true=x_eq_init+x_w; % init true sys state
x_ref_true=x_eq_ref+x_w;
tic;
x = x_w+x_eq_init; % real system input
for k = 1:(iterations)      
    fprintf('iteration no. %d/%d \n',k,iterations);
    
    
    
    % Implement first optimal control move and update plant states.
    [x, u] = getTransitionsTrue(x,x_w,r0,Kstable);
    
    % shift the output so that it's from the working point perspective
    % setpoint being [0;0;0;0]
    his = [x-x_w; u-r0]; % c: decision var; u-r0: delta u;
    
    % Save plant states for display.
    sysHistory = [sysHistory his]; 
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

figure;
plot_refs=plot(0:iterations, sysHistory(1:4,:), 'Linewidth',1.5);
grid on;
xlabel('iterations');
ylabel('responses');
plot_refs(1).Color = 'Yellow';
plot_refs(2).Color = 'Blue';
plot_refs(3).Color = 'Red';
plot_refs(4).Color = 'Green';

figure;
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1,'Marker','.');
grid on
xlabel('x1');
ylabel('x2');
title('State space');

