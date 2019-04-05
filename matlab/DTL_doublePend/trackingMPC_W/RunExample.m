%% TRACKING piecewise constant REFERENCE MPC example
clearvars;
close all;

%% Parameters
% Set the prediction horizon:
N = 2;

% The initial conditions
x = [0;...
     0];
%setpoint
xs = [0;...
      0.0];
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
% Simulation length (iterations)
iterations = 120;


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
K = -dlqr(A, B, Q, 1*R);
umin=-0.3;
max_admissible_ctrl_weight=1/(umin^2);
% K_t = -dlqr(A, B, Q, max_admissible_ctrl_weight*R);
K_t = -dlqr(A, B, Q, 10*R);
% K_t = -dlqr(A, B, 1000*Q, R); % empty U_bar set not good
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
% terminal steady state cost
T = 10*P;
%%
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================
u_min =[-0.3;-0.3]; u_max= [0.3;0.3];
x_min=[-5;-5]; x_max=[5;5];
w_min=[-0.1;-0.1]; w_max=[0.1;0.1]; % disturbance


F_u = [eye(m); -eye(m)]; h_u = [u_max; -u_min];
F_x = [eye(n); -eye(n)]; h_x = [x_max; -x_min];
F_g = [eye(n); -eye(n)]; h_g = [w_max; -w_min]; % uncertainty polytope
Xc=Polyhedron(F_x,h_x);
Uc=Polyhedron(F_u,h_u);
W=Polyhedron(F_g,h_g);

length_Fu = length(h_u);
length_Fx = length(h_x);
%%
%==========================================================================
% Compute maximally invariant set for tracking
%==========================================================================

disp('Computing and simplifying terminal set...');
L = (PSI - K*LAMBDA);
L0 = (PSI_0 - K*LAMBDA_0); % when being under inital disturbance
% contract the h values for the artificial steady state by a scalar λ ∈ (0, 1)
lambda=0.99;
% CONVEX POLYHAEDRON Wλ = {w = (x, θ) : (x, Kx + Lθ) ∈ Z, Mθθ ∈ λZ}.
F_w = [F_x zeros(length_Fx, m);
    zeros(length_Fx, n) F_x*LAMBDA; ...
    F_u*K, F_u*L; ...
    zeros(length_Fu, n) F_u*PSI];
h_w = [...
    h_x; ...
    lambda*(h_x - F_x*LAMBDA_0); ...
    h_u-F_u*L0; ...
    lambda*(h_u- F_u*PSI_0)];
F_w_N0 = F_w; h_w_N0 = h_w;
X_ext=Polyhedron(F_w,h_w);

% Compute usual maximally invariant set
syss = LTISystem('A', A+B*K);
poly = Polyhedron([K; -K; eye(2); -eye(2)], [u_max; -u_min; x_max; -x_min]);
iset1 = syss.invariantSet('X', poly);
MAIS_old=iset1.projection(1:n);
disp('Terminal set Polyhedron:');

% Compute new extended state maximally invariant set
Ak=[A+B*K B*L; zeros(m,n) eye(m)];
MPIS=compute_MPIS(X_ext,Ak);
MAIS=projection(MPIS,1:n); % Maximal Admissible Invariant set projected on X
F_w_N = MPIS.A; % Inequality description { x | H*[x; -1] <= 0 }
h_w_N = MPIS.b; % Inequality description { x | A*x <= b }

% Region of Attraction old
Xf0=MAIS_old;
XN0=ROA(params,Xf0,Xc,Uc,N);
% Region of attraction extended
Xf=X_ext.projection(1:2);
XN=ROA(params,Xf,Xc,Uc,N);

%% mRPI calculation
W=Polyhedron(W);
% S=Polyhedron(W.V);
% S.minVRep()
tic;
Z=reach_set(A+B*K_t,W,16);
toc
% figure;
% plot(Z);
F_z=Z.A; h_z=Z.b;
%% Robust running constraints
% count the length of the constraints on input, states, and uncertainty:
X_robust=Xc-Z;
figure;
plot([Xc,X_robust, Z]);
U_robust=Uc-K_t*Z;
figure;
plot([Uc,U_robust,K_t*Z]);

X_robust.minHRep(); U_robust.minHRep(); % simplify constraints
F_x=X_robust.A; h_x=X_robust.b;
F_u=U_robust.A; h_u=U_robust.b;

length_Fu = length(h_u);
length_Fx = length(h_x);

run_F = [F_x zeros(length_Fx,m);...
        zeros(length_Fu,n) F_u];
run_h = [h_x;h_u];
%% 
% extended Z constraints (for param theta) (Z x {0})
F_z_ext=[F_z zeros(length(h_z),m);
        zeros(m,n+m)];
h_z_ext = [h_z;
            zeros(m,1)];
Z_ext=Polyhedron(F_z_ext,h_z_ext);
Z_ext.minHRep();
% MAIS=projection(MPIS,1:n); % Maximal Admissible Invariant set projected on X

% x-theta constraints:
% F_xTheta = F_w_N;
% f_xTheta = h_w_N;
%% Robust terminal constraints
Xterm=MPIS;
MAIS_old=projection(Xterm,1:n); % Maximal Admissible Invariant set projected on X

XNr=ROA(params,MAIS_old,X_robust,U_robust,N);


F_xTheta = Xterm.A;
f_xTheta = Xterm.b;
% F_x = F_w_N(:, 1:n);
% F_theta = F_w_N(:,n+1:n+m);


%% Cost Calculation
% Start simulation
sysHistory = [x;u0(1:m,1)];
art_refHistory = Mtheta*theta0;
true_refHistory = xs;

for k = 1:(iterations)
    xs = set_ref(k);
    COSTFUN = @(var) costFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,xs,N,reshape(var(1:m),m,1),P,T,K,LAMBDA,PSI,params);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,N,K,LAMBDA,PSI,run_F,run_h,F_xTheta,f_xTheta,params);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    u = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    % Implement first optimal control move and update plant states.
    % Add disturbance
    x= getTransitions(x, u, params) + disturb(0.5*u_max,0.5*u_min); 
    his=[x; u];
    % Save plant states for display.
    sysHistory = [sysHistory his]; 
    art_refHistory = [art_refHistory art_ref];
    true_refHistory = [true_refHistory xs];
end


%% Plot
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
figure;
plot_refs=plot(0:iterations,art_refHistory(1,:), 0:iterations, true_refHistory(1,:),0:iterations,sysHistory(1,:),'Linewidth',2);
grid on;
legend({'art_{ref}','real_{ref}','x_1 response'},'Location','northeast'); 

xlabel('iterations');
% ylabel('references');
title('Artificial vs true reference vs state response');

plot_refs(1).LineStyle='--';
plot_refs(2).LineStyle='-.';
plot_refs(1).Color='green';
plot_refs(2).Color='red';
plot_refs(3).Color='b';

% set(gcf,'PaperPositionMode','auto')
% print('res','-dsvg','-r300') % set dpi to 300 and save in SVG
%%
figure;
% plot the sets
MAIS_old.plot('wire',1,'linewidth',2.5,'linestyle',':','color', 'blue'); 
hold on;
MAIS.plot('wire',1,'linewidth',2.5,'linestyle','-.','color', 'red'); 
hold on;
plot(XN,'wire',1,'linewidth',2.5,'linestyle','-','color', 'lightblue'); % ROA ext
hold on;
plot(XN0,'wire',1,'linewidth',2.5,'linestyle','--','color', 'green'); % ROA old
legend({'O_{\infty}(0)','X_f','X_N','X_N(O_{\infty}(0))'},'Location','northeast'); 
grid on;
xlabel('x_1');
ylabel('x_2');
title('Relevant sets');
% set(gcf,'PaperPositionMode','auto')
% print('sets','-dsvg','-r300') % set dpi to 300 and save in SVG
%%
figure;
Xfr=projection(Xterm,1:2);

XN.plot('wire',1,'linewidth',2.5,'linestyle','-'); hold on;
XNr.plot('wire',1,'linewidth',2.5,'linestyle','-.'); hold on;
Xf.plot('wire',1,'linewidth',2.5,'linestyle',':'); hold on;
Xfr.plot('wire',1,'linewidth',2.5,'linestyle','--');  % from L->R: bigger -> smaller set to have everything visible 
% plot the system state-space
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1.5,'Marker','o','color','k',  'MarkerSize',5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]); 
legend({'X_N','X_{N_{robust}}','X_f','X_{f_{robust}}','state'},'Location','northeast'); 
grid on;
xlabel('x_1');
ylabel('x_2');
title('State trajectory');
% set(gcf,'PaperPositionMode','auto')
% print('sstraj','-dsvg','-r300') % set dpi to 300 and save in SVG

%% Helper functions

%set reference depending on the iteration

function [xs] = set_ref(ct)
    if ct <=50
        xs=[3;0];
    elseif ct > 50 && ct <= 100
        xs=[-3;0];
    else
        xs=[0;0];
    end
end

function [w] = disturb(w_max,w_min)
    w = rand(2, 1).*(w_max - w_min)+w_min;
end

% function [w] = switching_diturb(w_max,w_min,ct)
%     if ct <=12
%         xs=[4;0];
%     elseif ct > 12 && ct <= 24
%         xs=[-4;0];
%     elseif ct > 12*2 && ct <= 24*3
%         xs=[-4;0];    
%     elseif ct > 12*3 && ct <= 24*4
%         xs=[-4;0];    
%     elseif ct > 12*5 && ct <= 24*6
%         xs=[-4;0];
%     end
% end