%% MPC for tracking piecewise constant references 
clearvars;
close all;

%% Parameters initialization
% Set the prediction horizon:
N = 3;
% Simulation length (iterations)
iterations = 100;

% The initial conditions
x = [0;...
     -2];
%setpoint
xs = [4.95;...
      0.0];
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');

% A non-square DLTI system - sampled double integrator
A = [1 1; 0 1];
B = [0.0 0.5; 1.0 0.5];
C = [1 0];
n = size(A,1);
m = size(B,2);
o = size(C,1);
Q = eye(n);
R = eye(m);

params.A=A;
params.B=B;
params.C=C;
params.n=n;
params.m=m;
params.Q=Q;
params.R=R;

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];

% Steady state parametrization matrices
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);

% For an initial disturbance guess
d_0 = [0,0]';
% Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + V*theta, theta in R^m
V_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = V_0(1:n);
PSI_0 = V_0(n+1:n+m);

%% Define a nominal feedback policy K and corresponding terminal cost
K = -dlqr(A, B, Q, R); % 'baseline' stabilizing feedback gain
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
% terminal steady state cost
T = 100*P;
%% Define polytopic constraints in H-representation
u_min =[-0.3;-0.3]; u_max= [0.3;0.3];
x_min=[-5;-5]; x_max=[5;5];
F_u = [eye(m); -eye(m)]; h_u = [u_max; -u_min];
F_x = [eye(n); -eye(n)]; h_x = [x_max; -x_min];
Xc=Polyhedron(F_x,h_x);
Uc=Polyhedron(F_u,h_u);

% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);

run_F = [F_x zeros(length_Fx,m);...
        zeros(length_Fu,n) F_u];
run_h = [h_x;h_u];
%% Compute maximally invariant set for tracking
L = (PSI - K*LAMBDA);
L0 = (PSI_0 - K*LAMBDA_0); % when being under inital disturbance
% contract the h values for the artificial steady state by a scalar \lambda
% to ensure that the invariant set is 
lambda=0.99;

% Convex polyhedrons of extended state for the DLTI system
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

% Compute the usual maximal admissible invariant set (MAIS)
syss = LTISystem('A', A+B*K);
poly = Polyhedron([K; -K; eye(2); -eye(2)], [u_max; -u_min; x_max; -x_min]);
iset1 = syss.invariantSet('X', poly);
MAIS_old=iset1.projection(1:n);
disp('Terminal set Polyhedron:');

% Compute new extended state maximally invariant set
Ak=[A+B*K B*L; zeros(m,n) eye(m)];
term_poly=compute_MPIS(X_ext,Ak);
MAIS=projection(term_poly,1:n); % Maximal Admissible Invariant set projected on X
F_w_N = term_poly.A; % Inequality description { x | H*[x; -1] <= 0 }
h_w_N = term_poly.b; % Inequality description { x | A*x <= b }

% Region of Attraction old
Xf0=MAIS_old;
XN0=ROA(params,Xf0,Xc,Uc,N);
iset2=syss.reachableSet('X',poly,'N',N,'direction','forward');

% Region of Attraction extended
Xf=X_ext.projection(1:2);
XN=ROA(params,Xf,Xc,Uc,N);

% x-theta constraints:
F_xTheta = F_w_N;
F_x = F_w_N(:, 1:n);
F_theta = F_w_N(:,n+1:n+m);
f_xTheta = h_w_N;

%% Cost Calculation
% Start simulation
sysHistory = [x;u0(1:m,1)];
art_refHistory = LAMBDA*theta0;
true_refHistory = xs;

for k = 1:(iterations)
    fprintf("Iterations : %d\n",k);
    xs = set_ref(k);
    COSTFUN = @(var) costFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,xs,N,reshape(var(1:m),m,1),P,T,K,LAMBDA,PSI,params);
    CONSFUN = @(var) constraintsFunction(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),x,N,K,LAMBDA,PSI,run_F,run_h,F_xTheta,f_xTheta,params);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    u = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    % Implement first optimal control move and update plant states.
    x= getTransitions(x, u, params); %-K*(x-xs)+u_opt
    his=[x; u];
    % Save plant states for display.
    sysHistory = [sysHistory his]; %#ok<*AGROW>
    art_refHistory = [art_refHistory art_ref(1:n)];
    true_refHistory = [true_refHistory xs];
end

%% Plot

disp('Generating plots...');
figure;
for i=1:n+m
        subplot(n+m,1,i);
        plot(0:iterations,sysHistory(i,:),'Linewidth',1); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = [' x_{',num2str(i),'}'];
        else
            lableTex = [' u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
 end
suptitle('System response');
% 
figure;
plot_refs=plot(0:iterations,art_refHistory(1,:), 0:iterations, true_refHistory(1,:),0:iterations,sysHistory(1,:),'Linewidth',1.5);
grid on;
legend({'art_{ref}','real_{ref}','x_1 response'},'Location','southeast'); 
xlabel('iterations');
title('x_1 response');

plot_refs(1).LineStyle='--';
plot_refs(2).LineStyle='-.';
plot_refs(1).Color='green';
plot_refs(2).Color='red';
plot_refs(3).Color='b';

%%
figure;
% plot the sets
MAIS_old.plot('wire',1,'linewidth',2.5,'linestyle',':','color', 'blue'); 
hold on;
MAIS.plot('wire',1,'linewidth',2.5,'linestyle','-.','color', 'red'); 
hold on;
plot(XN,'wire',1,'linewidth',2.5,'linestyle','-','color', 'lightblue'); % ROA ext
hold on;
plot(XN0,'wire',1,'linewidth',2.5,'linestyle','-','color', 'green'); % ROA old
legend({'O_{\infty}(0)','X_f','X_N','X_N(O_{\infty}(0))'},'Location','southwest'); 
grid on;
xlabel('x_1');
ylabel('x_2');
title('Relevant sets');

%%
figure;
plot([XN,MAIS]); % from L->R: bigger->smaller set to have everything visible 
hold on;
% plot the system state-space
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1.5,'Marker','o','color','k',  'MarkerSize',5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]); 
legend({'X_N','X_f','state'},'Location','southwest'); 
grid on;
xlabel('x_1');
ylabel('x_2');
title('State trajectory');

%% Helper functions

%set reference depending on the iteration
function [xs] = set_ref(ct)
    if ct <=30 
        xs=[4.95;0];
    elseif ct > 30 && ct <= 60
        xs=[-5.5;0];
    elseif ct > 60 && ct <= 90
        xs=[2;0];
    else
        xs=[0;0];
    end
end