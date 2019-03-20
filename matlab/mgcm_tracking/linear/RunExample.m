%% TRACKING piecewise constant REFERENCE MPC example
close all;
clearvars;

%% Parameters
% Set the prediction horizon:
N = 10;
% Simulation length (iterations)
iterations = 60;

% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;
%% Closed-Loop Simulation
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
% Mtheta = [1, 0, 0, 0; 0, 1, 1, -2]';
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);
% 
Q = eye(n);
R = eye(m);
[P, e, K] = dare(A,B,Q,R);
K=[-3.69438241357512	-1.78296525690465	-0.271009055947981	0.00595919621278219];
T = 100*P; 

%% Calculating constraints
% Optimization variable bounds (usually control) not constraints per say
LB = [ones(m*N,1)*u_min; ones(m,1)*(-inf)];
UB = [ones(m*N,1)*u_max; ones(m,1)*(+inf)];
u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];

ALPHA = 0.99; 

Mtheta =[LAMBDA' PSI']';

L = [K eye(m)]*Mtheta;

sysStruct.A=[A-B*K ,     B*L;...
            zeros(m,n), eye(m)];
sysStruct.B=zeros(n+m,m);
sysStruct.C=zeros(m,n+m);
sysStruct.D=zeros(m,m);

% Q = diag([1,1]);
% R = diag([1,1]);
% [P, e, K] = dare(A,B,Q,R);
umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max;throttle_max;throttle_rate_max]; 
xmin = [mflow_min; prise_min;throttle_min;throttle_rate_min];
run_F = [eye(m+n);eye(m+n)];
run_h= [xmax;umax;-xmin;-umin];

%% terminal constraints from invariant set calculation
% sysStruct.xmax = [xmax; inf(m,1)]*ALPHA;
% sysStruct.xmin = [xmin;-inf(m,1)]*ALPHA;
% sysStruct.umax = umax*ALPHA;
% sysStruct.umin= umin*ALPHA;
% % sysStruct.umax = umax*ALPHA;
% % sysStruct.umin= umin*ALPHA;
% % sysStruct.x.penalty.weight=Q;
% % sysStruct.u.penalty.weight=R;
% 
% system = LTISystem(sysStruct);
% InvSet2 = system.invariantSet(); % InvSet2 is a polyhaeder
% % extracting H-representation
% term_F=InvSet2.A;
% term_h=InvSet2.b;
% % InvSet2.plot();
% 
% % project the 4D case to a 2D one
% % MAI=projection(InvSet2,1:2); % Maximal Admissible Invariant set
% % plot(MAI);
term_F = []; term_h = [];
%% Start simulation
sysHistory = [x;u0(1:m,1)];
% art_refHistory = LAMBDA*theta0;
art_refHistory =  0;
true_refHistory = xs;

for k = 1:(iterations)
    COSTFUN = @(var) costFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,xs,N,reshape(var(1:m),m,1),Q,R,P,T,K,LAMBDA,PSI);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,N,K,LAMBDA,PSI,run_F,run_h,term_F,term_h);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],LB,UB,CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
%     c=0;
    art_ref = Mtheta*theta_opt;
    % Implement first optimal control move and update plant states.
    x= getTransitions(x, c, K, theta_opt, LAMBDA, PSI); %-K*(x-xs)+u_opt
    his=[x; c]; % c = delta u
    % Save plant states for display.
    sysHistory = [sysHistory his]; 
    art_refHistory = [art_refHistory art_ref(1:m)];
    true_refHistory = [true_refHistory xs];
end

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

