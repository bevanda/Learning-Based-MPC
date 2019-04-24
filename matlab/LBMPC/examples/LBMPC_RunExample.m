clearvars;  
addpath('../'); 
%% INIT CONTROLLER DESIGN
syms u ... % control input
    x1 ... % mass flow
    x2 ... % pressure rise
    x3 ... % throttle opeining
    x4 ... % throttle opening rate
    ;
wn=sqrt(1000); % resonant frequency
zeta=1/sqrt(2); % damping coefficient
beta=1; % constant >0
x2_c=0; % pressure constant

%% Continous time state-space model of the Moore-Greitzer compressor model

f1 = -x2+x2_c+1+3*(x1/2)-(x1^3/2); % mass flow rate
f2 = (x1+1-x3*sqrt(x2))/(beta^2); % pressure rise rate
f3 = x4; % throttle opening rate
f4 = -wn^2*x3-2*zeta*wn*x4+wn^2*u; % throttle opening acceleration

%% Linearisation around the equilibrium [0.5 1.6875 1.1547 0]'

A = jacobian([f1,f2, f3, f4], [x1, x2, x3, x4]);
B = jacobian([f1,f2, f3, f4], [u]);

% equilibrium params
x1 = 0.5;
x2 = 1.6875;
x3 = 1.1547;
x4 = 0;
equili = [x1 x2 x3 x4];
init_cond = [x1-0.35, x2-0.4, x3, 0];  % init condition
% print the matrices in the cmd line
A = eval(A);
B = eval(B);
% C = [A(1, :); A(2,:)] % choose f1 and f2 as outputs 
C = eye(4);
D = zeros(4,1);
n=size(A,2);

% Visualise the poles and zeros of the continuous system
[b,a]=ss2tf(A,B,C,D);
sys = tf([b(1,:)],[a]);
% sys2 = tf([b(2,:)],[a]);
% figure;
% pzmap(sys);
% grid on;
% pzmap(sys2);
%% Exact discretisation

dT = 0.01; % sampling time

Ad = expm(A*dT);
Bd = (Ad-eye(n))*inv(A)*B;
Cd=C;
Dd=D;
Td=dT;
Ts=dT;
e = eig(Ad);
% figure;
% sys = idss(Ad,Bd,Cd,Dd,'Ts',dT);
% pzmap(sys);

%% System stabilisation /w feedback matrix K to place poles near Re(p_old) inside unit circle
switch dT
    case 0.05
        p=[0.13, 0.16, 0.9, 0.95]; % dT=0.05
    case 0.02
        p=[0.55, 0.6, 0.9, 0.95]; % dT=0.02
    case 0.01
        p=[0.75, 0.78, 0.98, 0.99]; % dT=0.01
end

[K,prec,message] = place(Ad,Bd,p) %nominal feedback matrix
Kstabil=-K;

AK = Ad+Bd*Kstabil;
e = eig(AK);

Q = eye(4); R=1;
P = dare(AK,Bd,Q,R);
Klqr= -dlqr(Ad,Bd,Q,R);
% figure;
% sys = idss(AK,zeros(4,1),Cd,Dd,'Ts',dT);
% pzmap(sys);



%% Parameters
% Horizon length
N=10;
% Simulation length (iterations)
iterations = 10/dT;

%% Discrete time nominal model of the non-square LTI system for tracking
A = Ad;
B = Bd;
C = Cd;
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

Klqr = -dlqr(A, B, Q, R);
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
% state_uncert = [0.01;0.01;0.01;0.01]; % test
%%
%  Shift the constraints for the linearised model for the value of the
%  working point
x_w = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
u_w = x_w(3);

%% testing /w oracle
% shrnik=0.05;
shrnik=0.0;
% Shift the abs system constraints w.r.t. to the linearisation point
F_u = [eye(m); -eye(m)]; h_u = [umax-u_w-shrnik; -umin+u_w+shrnik];
F_x = [eye(n); -eye(n)]; h_x = [xmax-x_w-shrnik; -xmin+x_w+shrnik];
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
    [F_x_d, h_x_d] = double(temp);
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
maxadm_controlweight = 10; % r_i as the inverse of the square of the maximum permissible value for the corresponding u_i
K_t = -dlqr(A, B, Q, maxadm_controlweight*R);
lambda=0.99; % λ ∈ (0, 1), λ can be chosen arbitrarily close to 1, the obtained
% invariant set can be used as a reliable polyhedral approximation to the maximal invariant set 
disp('Computing and simplifying terminal set...');
% extended state constraints
L=(PSI - K_t*LAMBDA);
L0=(PSI_0 - K_t*LAMBDA_0);
F_w = [ F_x zeros(length_Fx, m);
        zeros(length_Fx, n) F_x*LAMBDA; ...
        F_u*K_t, F_u*L; ...
        zeros(length_Fu, n) F_u*PSI; ...
        F_x_d*(A+B*K_t) F_x_d*B*L];
h_w = [ h_x; ...
        lambda*(h_x - F_x*LAMBDA_0); ...
        h_u - F_u*L0; ...
        lambda*(h_u - F_u*PSI_0); ...
        h_x_d - F_x_d*B*(PSI_0-K_t*LAMBDA_0)];

   

% disturbance constraints of the extended state 
F_g_w = [F_g zeros(length_Fg,m); ...
        zeros(m, n) eye(m); ...
        zeros(m, n) -eye(m)];
h_g_w = [h_g; ...
        zeros(2*m,1)];
    

% calculating the robust positively invariant set    
[F_w_N0, h_w_N0] = pdiff(F_w, h_w, F_g_w, h_g_w);


% Simplify the constraints
term_poly = polytope(F_w_N0, h_w_N0);
[F_w_N, h_w_N] = double(term_poly);
% term_poly = Polyhedron(F_w_N01, h_w_N01);
% term_poly.minHRep(); % simplifying the polyhedron/constraints
% F_w_N1 = term_poly.A; h_w_N1 = term_poly.b;
disp('Terminal set Polyhedron:');
term_poly
% term_poly2=polytope(F_w,h_w)

%% Start simulation

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];

sysHistory = [x_eq_init;u0(1:m,1)];

art_refHistory =  0;
true_refHistory = x_eq_ref;
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
x = x_w+x_eq_init; % true sistem init state
% init states from models used in MPC
xl=x_eq_init;
xo=x_eq_init;
% init data from estimation
data.X=zeros(3,1);
data.Y=zeros(4,1);
% sysHistoryO=[x_eq_init; 0];
%% Run LBMPC
tic;
for k = 1:(iterations)      
    fprintf('iteration no. %d/%d \n',k,iterations);
    if k>1
        % DATA ACQUISTION 
        X=[x(1:2)-x_w(1:2); u-u_w]; %[δphi;δpsi;δu]
        switch dT
            case 0.01
                Y=((x_k1-x_w)-(A*(x-x_w)+B*(u-u_w))); %[δx_true-δx_nominal]
            otherwise
                Y=-((x_k1-x_w)-(A*(x-x_w)+B*(u-u_w))); %[δx_nominal-δx_true]
        end
        % update state vars for estimation
        x=x_k1;
%         xl=xn_k1;
        % get iterations
        q=200; % moving window of q datapoints 
        data=update_data(X,Y,q,k,data);
        
        % get the real state w.r.t. equilibrium
        x_eq=x-x_w;
    else
        x_eq=x_eq_init;
    end
    
    % SOLVE THE OPTIMAL CONTROL PROBLEM
    COSTFUN = @(var) costLBMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        x_eq,x_eq_ref,N,reshape(var(1:m),m,1),Q,R,P,T,Kstabil,x_w,u_w,LAMBDA,PSI,data,dT);
    CONSFUN = @(var) constraintsLBMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        x_eq,N,Kstabil,F_x,h_x,F_u,h_u,F_w_N,h_w_N,F_x_d,h_x_d);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    % Apply control to system and models
    % Implement first optimal control move and update plant states.
    [x_k1, u] = transitionTrue(x,c,x_w,u_w,Kstabil,dT); % plant   
    
    % Save state data for plotting w.r.t. work point x_w
    % shift the output so that it's from the working point perspective
    % setpoint being [0;0;0;0]
    [x_wp, u_wp]=wp_shift(x,x_w,u,u_w);
    his = [x_wp; u_wp]; 
    % Save plant states for display.
    sysHistory = [sysHistory his]; %#ok<*AGROW>
%     sysHistoryO = [sysHistoryO hisO]; %#ok<*AGROW>
    art_refHistory = [art_refHistory art_ref(1:m)];
    true_refHistory = [true_refHistory x_eq_ref];
    
end
toc

%% PLOT
figure;
subplot(n+m,1,1);
plot(0:iterations,sysHistory(1,:),'Linewidth',1.5); hold on;
% plot(0:iterations,sysHistoryO(1,:),'Linewidth',1.5); 
grid on
xlabel('iterations');
ylabel('x1');
title('mass flow');
subplot(n+m,1,2);
plot(0:iterations,sysHistory(2,:),'Linewidth',1.5); hold on;
% plot(0:iterations,sysHistoryO(2,:),'Linewidth',1.5); 
grid on
xlabel('iterations');
ylabel('x2');
title('pressure rise');
subplot(n+m,1,3);
plot(0:iterations,sysHistory(3,:),'Linewidth',1.5); hold on;
% plot(0:iterations,sysHistoryO(3,:),'Linewidth',1.5); 
grid on
xlabel('iterations');
ylabel('x3');
title('throttle');
subplot(n+m,1,4);
plot(0:iterations,sysHistory(4,:),'Linewidth',1.5); hold on;
% plot(0:iterations,sysHistoryO(4,:),'Linewidth',1.5); 
grid on
xlabel('iterations');
ylabel('x4');
title('throttle rate');
subplot(n+m,1,5);
plot(0:iterations,sysHistory(5,:),'Linewidth',1.5); hold on;
% plot(0:iterations,sysHistoryO(5,:),'Linewidth',1.5); 
grid on
xlabel('iterations');
ylabel('u');
title('Sys input');
legend({'System response'});


figure;
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1.5,'LineStyle','-'); hold on;
grid on
xlabel('x1');
ylabel('x2');
title('State space');