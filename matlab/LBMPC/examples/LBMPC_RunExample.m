clearvars;  
addpath('../'); 
addpath('../misc/'); 
addpath('../models/'); 

%% Parameter initialization
% Generate the DLTI verison of the continouous Moore-Greitzer Compressor 
% model (MGCM)
[A,B,C,D,Ts]=mgcmDLTI();
n = size(A,1); % num of states
m = size(B,2); % num of inputs
o = size(C,1); % num of outputs

% Obtaining all needed matrices for the optimal control problem (OCP)
[Kstabil,Klqr,Q,R,P,T,Mtheta,LAMBDA,PSI,LAMBDA_0,PSI_0]=matOCP(A,B,C,n,m,o);

%% OCP setting
N=60; % Horizon length
iterations = 10/Ts; % Simulation length (iterations)

% The initial conditions w.r.t. to the linearisation/working point
x_wp_init = [-0.35;...
    -0.4;...
    0.0;...
    0.0];
% setpoint
x_wp_ref = [0.0;...
      0.0;...
      0.0;...
      0.0];
 
% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;

umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];

% Uncertainty bound
state_uncert = [0.02;5e-04;0;0]; % obtained form the maximal 
% linearization error Lagrange Error Bound + estimation error toleranace

% Working point (wp)
x_wp = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
u_wp = x_wp(3);

[F_x,h_x, F_u,h_u, F_w_N,h_w_N, F_x_d,h_x_d]...
    =getCONSPOLY(xmax,xmin,umax,umin,state_uncert,x_wp,u_wp,m,n,...
    A,B,Q,R,LAMBDA,PSI,LAMBDA_0,PSI_0);
%% Simulation setup

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];

sysHistory = [x_wp_init;u0(1:m,1)];

art_refHistory =  0;
true_refHistory = x_wp_ref;
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
x = x_wp+x_wp_init; % true sistem init state
% init states from models used in MPC
xl=x_wp_init;
xo=x_wp_init;
% init data from estimation
data.X=zeros(3,1);
data.Y=zeros(4,1);

%% Run LBMPC
tic;
for k = 1:(iterations)      
    fprintf('iteration no. %d/%d \n',k,iterations);
    if k>1
        % DATA ACQUISTION 
        X=[x(1:2)-x_wp(1:2); u-u_wp]; %[δphi;δpsi;δu]
        switch Ts
            case 0.01
                Y=((x_k1-x_wp)-(A*(x-x_wp)+B*(u-u_wp))); %[δx_true-δx_nominal]
            otherwise
                Y=-((x_k1-x_wp)-(A*(x-x_wp)+B*(u-u_wp))); %[δx_nominal-δx_true]
        end
        % update state vars for estimation
        x=x_k1;
%         xl=xn_k1;
        % get iterations
        q=100; % moving window of q datapoints 
        data=update_data(X,Y,q,k,data);
        
        % get the real state w.r.t. wpuilibrium
        dx=x-x_wp;
    else
        dx=x_wp_init;
    end
    
    % SOLVE THE OPTIMAL CONTROL PROBLEM
    COSTFUN = @(var) costLBMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        dx,x_wp_ref,N,reshape(var(1:m),m,1),Q,R,P,T,Kstabil,x_wp,u_wp,LAMBDA,PSI,data,Ts);
    CONSFUN = @(var) constraintsLBMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        dx,N,Kstabil,F_x,h_x,F_u,h_u,F_w_N,h_w_N,F_x_d,h_x_d);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;
    % Apply control to system and models
    % Implement first optimal control move and update plant states.
    [x_k1, u] = transitionTrue(x,c,x_wp,u_wp,Kstabil,Ts); % plant   
    
    % Save state data for plotting w.r.t. work point x_w
    % shift the output so that it's from the working point perspective
    % setpoint being [0;0;0;0]
    [dx, du]=wp_shift(x,x_wp,u,u_wp);
    his = [dx; du]; 
    % Save plant states for display.
    sysHistory = [sysHistory his]; %#ok<*AGROW>
    art_refHistory = [art_refHistory art_ref(1:m)];
    true_refHistory = [true_refHistory x_wp_ref];
    
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