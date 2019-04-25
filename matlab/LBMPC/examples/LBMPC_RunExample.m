clearvars;  
addpath('../'); 
% Generate the DLTI verison of the Moore-Greitzer Compressor model (mgcm)
[A,B,C,D,Ts]=mgcmDLTI();
n = size(A,1); % num of states
m = size(B,2); % num of inputs
o = size(C,1); % num of outputs

% System stabilisation /w feedback matrix K and
% state space parametrization
[Kstabil,Klqr,Q,R,P,T,Mtheta,LAMBDA,PSI,LAMBDA_0,PSI_0]=matOCP(A,B,C,n,m,o);


%% Parameters
% Horizon length
N=10;
% Simulation length (iterations)
iterations = 10/Ts;
%% Simulation setup

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

state_uncert = [0.02;5e-04;0;0]; % from max lin error from TaylorRemainder (Lagrange Error Bound)
%%
%  Shift the constraints for the linearised model for the value of the
%  working point
x_w = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
u_w = x_w(3);
%%
% Shift the abs system constraints w.r.t. to the linearisation point
F_u = [eye(m); -eye(m)]; h_u = [umax-u_w; -umin+u_w];
F_x = [eye(n); -eye(n)]; h_x = [xmax-x_w; -xmin+x_w];
F_g = [eye(n); -eye(n)]; h_g = [state_uncert; state_uncert]; % uncertainty polytope
% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);
length_Fg = length(h_g);
%% State constraints
temp = polytope(F_x, h_x) - polytope(F_g, h_g);
[F_x_d, h_x_d] = double(temp);
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
        switch Ts
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
        x_eq,x_eq_ref,N,reshape(var(1:m),m,1),Q,R,P,T,Kstabil,x_w,u_w,LAMBDA,PSI,data,Ts);
    CONSFUN = @(var) constraintsLBMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        x_eq,N,Kstabil,F_x,h_x,F_u,h_u,F_w_N,h_w_N,F_x_d,h_x_d);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    % Apply control to system and models
    % Implement first optimal control move and update plant states.
    [x_k1, u] = transitionTrue(x,c,x_w,u_w,Kstabil,Ts); % plant   
    
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