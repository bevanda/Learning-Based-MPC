clearvars;
close all;
%% Robust MPC for tracking piecewise constant references (without disturbance cancelation)
% Set the prediction horizon:
N = 3;
% The initial conditions
x = [0;...
     0];
%setpoint
xs = [0;...
      0.0];
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
% Simulation length (iterations)
iterations = 150;


A = [1 1; 0 1];
B = [0.0 0.5; 1.0 0.5];
C = [1 0];
n = size(A,1);
m = size(B,2);
o = size(C,1);
% Optimization variable bounds not constraints per say

params.A=A;
params.B=B;
params.C=C;
params.n=n;
params.m=m;

Q = eye(n);
R = eye(m);
params.Q=Q;
params.R=R;
u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];

% MN = [Mtheta; 1, 0];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);


% init disturbance estimate
d_0 = [0,0]';
% Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + V*theta, theta in R^m
V_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = V_0(1:n);
PSI_0 = V_0(n+1:n+m);


%%
% Define a nominal feedback policy K and corresponding terminal cost
% 'baseline' stabilizing feedback law
R=10*R;
K = -dlqr(A, B, Q, R);
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
%% terminal steady state cost
T = 100*P;

% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g

u_min =[-0.3;-0.3]; u_max= [0.3;0.3];
x_min=[-5;-5]; x_max=[5;5];
w_min=[-0.1;-0.1]; w_max=[0.1;0.1]; % disturbance


F_u = [eye(m); -eye(m)]; h_u = [u_max; -u_min];
F_x = [eye(n); -eye(n)]; h_x = [x_max; -x_min];
F_d = [eye(n); -eye(n)]; h_d = [w_max; -w_min]; % uncertainty polytope
Xc=Polyhedron(F_x,h_x);
Uc=Polyhedron(F_u,h_u);
W=Polyhedron(F_d,h_d);

length_Fu = length(h_u);
length_Fx = length(h_x);


L = (PSI - K*LAMBDA);
L0 = (PSI_0 - K*LAMBDA_0); % when being under inital disturbance
% contract the h values for the artificial steady state by a scalar λ ∈ (0, 1)
lambda=0.95;

F_w = [F_x zeros(length_Fx, m);
    zeros(length_Fx, n) F_x*LAMBDA; ...
    F_u*K, F_u*L; ...
    zeros(length_Fu, n) F_u*PSI];
h_w = [...
    h_x; ...
    lambda*(h_x - F_x*LAMBDA_0); ...
    h_u-F_u*L0; ...
    lambda*(h_u- F_u*PSI_0)];
X_ext=Polyhedron(F_w,h_w);

% MPIS
Ak=[A+B*K B*L; ...
    zeros(m,n) eye(m)];
MPIS=compute_MPIS(X_ext,Ak);
% ROA 
Xf=MPIS.projection(1:n);
XN=ROA(params,Xf,Xc,Uc,N);

% mRPI set calculation;
tic;
Z=calc_mRPIS(A+B*K,F_d,h_d,5e-6);
toc
F_z=Z.A; h_z=Z.b;

X_c_r=Xc-Z;
U_c_r=Uc-K*Z;

X_c_r.minHRep(); U_c_r.minHRep(); % simplify constraints
F_x_z=X_c_r.A; h_x_z=X_c_r.b;
F_u_z=U_c_r.A; h_u_z=U_c_r.b;

length_Fuz = length(h_u_z);
length_Fxz = length(h_x_z);
run_F = [F_x_z zeros(length_Fxz,m);...
        zeros(length_Fuz,n) F_u_z];
run_h = [h_x_z;h_u_z];

% extended Z constraints (for param theta) (Z x {0})
Ftheta=[eye(m); -eye(m)];
length_Ftheta=length(Ftheta);
length_Fz=length(F_z);
F_z_ext=[      F_z   ,   zeros(length(h_z),m);
        zeros(length_Ftheta,m), Ftheta];
h_z_ext = [h_z;
            zeros(length_Ftheta,1)];
Z_ext=Polyhedron(F_z_ext,h_z_ext);
Z_ext.minHRep();

% MPIS_Z=MPIS-Z_ext;
% MPIS_Z.minHRep();
MPIS_Z=compute_MPIS(X_ext-Z_ext,Ak);
F_w_N=MPIS_Z.A; h_w_N=MPIS_Z.b;
% ROA robust
Xfr=MPIS_Z.projection(1:n);
XNr=ROA(params,Xfr,X_c_r,U_c_r,N);

% Start simulation
sysHistory = [x;u0(1:m,1)];
art_refHistory = Mtheta*theta0;
true_refHistory = xs;

for k = 1:(iterations)
    fprintf("Iterations : %d\n",k);
    xs = set_ref(k);
    COSTFUN = @(var) costFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,xs,N,reshape(var(1:m),m,1),P,T,K,LAMBDA,PSI,params);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,N,K,LAMBDA,PSI,run_F,run_h,F_w_N,h_w_N,params);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    u = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    % Implement first optimal control move and update plant states.
    % Add disturbance
    x= getTransitions(x, u, params) +... %disturb(w_max,w_min);
     1*switching_diturb(w_max(1),w_min(1),w_max(2),w_min(2),k);
    his=[x; u];
    % Save plant states for display.
    sysHistory = [sysHistory his]; %#ok<*AGROW> 
    art_refHistory = [art_refHistory art_ref]; 
    true_refHistory = [true_refHistory xs];
end


disp('Generating plots...');

figure;
subplot(2,1,1);   
plot(0:iterations,sysHistory(1:n,:),'Linewidth',2); hold on;
plot(0:iterations,art_refHistory(1:n,:),'Linewidth',1.5,'LineStyle','--');
grid on;
legend({'ref x_1','ref x_2','x_1 response','x_2 response'},'Location','northeast');
xlabel('iterations');
ylabel('x');
% title('states');
subplot(2,1,2);   
plot(0:iterations,sysHistory(n+1:n+m,:),'Linewidth',2); hold on;
plot(0:iterations,art_refHistory(n+1:n+m,:),'Linewidth',1.5,'LineStyle','--');
grid on;
legend({'ref u_1','ref u_2','u_1 response','u_2 response'},'Location','northeast');
xlabel('iterations');
ylabel('u');
    
%%
% x1 plot
figure;
plot_refs=plot(0:iterations,art_refHistory(1,:), 0:iterations, true_refHistory(1,:),0:iterations,sysHistory(1,:),'Linewidth',2);
grid on;
legend({'art_{ref}','real_{ref}','x_1 response'},'Location','northeast'); 

xlabel('iterations');
title('Artificial vs true reference vs state response');

plot_refs(1).LineStyle='--';
plot_refs(2).LineStyle='-.';
plot_refs(1).Color='green';
plot_refs(2).Color='red';
plot_refs(3).Color='b';

% set plots
figure;
plot(Z);
xlabel('x_1');
ylabel('x_2');
title('mRPI')
figure;
plot([Xc,X_c_r, Z]);
xlabel('x_1');
ylabel('x_2');
legend({'X','X_{robust}','Z'},'Location','northeast'); 
figure;
plot([Uc,U_c_r,K*Z]);
xlabel('u_1');
ylabel('u_2');
legend({'U','U_{robust}','KZ'},'Location','northeast');

figure;
Xf=projection(Xf,1:2);
Xfr=projection(Xfr,1:2);
XN=projection(XN,1:2); 
XNr=projection(XNr,1:2); 
XN.plot('wire',1,'linewidth',2.5,'linestyle','-'); hold on;
XNr.plot('wire',1,'linewidth',2.5,'linestyle','-.'); hold on;
Xf.plot('wire',1,'linewidth',2.5,'linestyle',':'); hold on;
Xfr.plot('wire',1,'linewidth',2.5,'linestyle','--');   % from L->R: bigger -> smaller set to have everything visible 
legend({'X_{N}','X_{N_{robust}}','X_f','X_{f_{robust}}'},'Location','northeast'); 
grid on;
xlabel('x_1');
ylabel('x_2');
title('Relevant sets');

%% Helper functions
function [xs] = set_ref(ct)
    if ct <=50
        xs=[4;0];
    elseif ct > 50 && ct <= 100
        xs=[-4;0];
    else
        xs=[0;0];
    end
end

function [w] = disturb(w_max,w_min)
    w = rand(2, 1).*(w_max - w_min)+w_min;
end

function [w] = switching_diturb(w1_max,w1_min,w2_max,w2_min,ct)
    if ct <=12
        w=[w1_max;w2_max];
    elseif ct > 12 && ct <= 24
        w=[w1_min;w2_max];
    elseif ct > 12*2 && ct <= 12*3
        w=[w1_max;w2_min];    
    elseif ct > 12*3 && ct <= 12*4
        w=[w1_max;w2_max];    
    elseif ct > 12*4 && ct <= 12*5
        w=[w1_min;w2_max];
    elseif ct > 12*5 && ct <= 12*6
        w=[w1_min;w2_min];
    elseif ct > 12*6 && ct <= 12*7
        w=[w1_max;w2_min];
    elseif ct > 12*7 && ct <= 12*8
        w=[w1_max;w2_max];    
    elseif ct > 12*8 && ct <= 12*9
        w=[w1_min;w2_max];    
    elseif ct > 12*49 && ct <= 12*10
        w=[w1_max;w2_max];
    elseif ct > 12*10 && ct <= 12*11
        w=[w1_min;w2_min];
    else %ct > 12*11 && ct <= 12*12
        w=[w1_max;w2_min];
    end
end