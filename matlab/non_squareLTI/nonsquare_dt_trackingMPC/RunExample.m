%% TRACKING piecewise constant REFERENCE MPC example
close all;
% clearvars;
%% Calculate the invariant set for tracking using MPT ..................

%% Parameters
% Set the prediction horizon:
N = 3;

%% Closed-Loop Simulation
% The initial conditions
x = [0;...
    -2];
%setpoint
xs = [4.95;...
      0.0];

options = optimoptions('fmincon','Algorithm','sqp','Display','none');

% Simulation time in num iter
iterations = 90;
% Optimization variable bounds (usually control) not constraints per say
u_min =-0.3; u_max= 0.3;
LB = [ones(2*N,1)*u_min; ones(2,1)*(-inf)];
UB = [ones(2*N,1)*u_max; ones(2,1)*(+inf)];
% LB =[];
% UB = [];

A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
% Count number of states n, number of inputs m, number of outputs o:
n = size(A,1);
m = size(B,2);
o = size(C,1);

Ntheta = [1, 0];
% MN = [Mtheta; Ntheta];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
% Mtheta = [1, 0, 0, 0; 0, 1, 1, -2]';
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);

Q = diag([1,1]);
R = diag([1,1]);
[P, e, K] = dare(A,B,Q,R);
T = 100*P;

% for inv set for tracking calculation
ALPHA = 0.99; 
L = [K eye(m)]*[LAMBDA' PSI']';
A_w = [A-B*K, B*L;...       % extended state w=(x,theta)
       zeros(n,m), eye(m)];

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];
%% Cost Calculation
% Start simulation
sysHistory = [x;u0(1:2,1)];
art_refHistory = LAMBDA*theta0;
true_refHistory = xs;
% J=costFunction(reshape(opt_var(1:end-2),2,N),reshape(opt_var(end-1:end),2,1),x,xs,N,reshape(opt_var(1:2),2,1),reshape(opt_var(end-1:end),2,1),P,T,LAMBDA,PSI);
for k = 1:(iterations)
    xs = set_ref(k);
    % opt_var must be a vector!
    COSTFUN = @(var) costFunction(reshape(var(1:end-2),2,N),reshape(var(end-1:end),2,1),x,xs,N,reshape(var(1:2),2,1),reshape(var(end-1:end),2,1),P,T,K,LAMBDA,PSI);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-2),2,N),reshape(var(end-1:end),2,1),x,N,K,xs);
%     [c, ceq]=CONSFUN(opt_var)
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],LB,UB,CONSFUN,options);    
    theta_opt = reshape(opt_var(end-1:end),2,1);
    u_opt = reshape(opt_var(1:2),2,1);
    art_ref = Mtheta*theta_opt;
    % Implement first optimal control move and update plant states.
    x = getTransitions(x, u_opt, K, xs); %-K*(x-xs)+u_opt
    his=[x; u_opt];
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
plot_refs=plot(0:iterations,art_refHistory(1,:), 0:iterations, true_refHistory(1,:),0:iterations,sysHistory(1,:),'Linewidth',1);
grid on
xlabel('iterations');
% ylabel('references');
title('Artificial vs true reference vs state response');
legend({'artifical reference','real reference', 'state response'},'Location','northeast')
plot_refs(1).Marker='.';
plot_refs(2).Marker='.';
plot_refs(1).Color='Green';
plot_refs(2).Color='Red';
plot_refs(3).Color='Blue';

figure;
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1,'Marker','.');
grid on
xlabel('x1');
ylabel('x2');
title('State space');

%% Help functions
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