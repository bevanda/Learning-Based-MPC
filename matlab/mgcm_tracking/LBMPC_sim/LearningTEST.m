%% TRACKING piecewise constant REFERENCE MPC example
close all;
clearvars;

%% Parameters
% Horizon length
N=20;
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
x_w = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
r0 = x_w(3);
Q = eye(n);
R = eye(m);
%%
%==========================================================================
% Define a nominal feedback policy K and corresponding terminal cost
% 'baseline' stabilizing feedback law
K = -dlqr(A, B, Q, R);


%%
% Kstable=[+3.0741 2.0957 0.1197 -0.0090]; % K stabilising gain from the papers
Kstable=-[3.0742   -2.0958   -0.1194    0.0089];

u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values
opt_var = [u0; theta0];
%% Start simulation
sysHistory = [x_eq_init;u0(1:m,1)];
sysHistoryL=sysHistory;
sysHistoryO=sysHistory;

art_refHistory =  0;
true_refHistory = x_eq_ref;
options = optimoptions('fmincon','Algorithm','sqp','Display','notify');
x_init_true=x_eq_init+x_w; % init true sys state
x_ref_true=x_eq_ref+x_w;
x = x_w+x_eq_init; % real system input
xl=x_eq_init;
xo=x_eq_init;
data.X=zeros(3,1);
data.Y=zeros(4,1);
tic;
for k = 1:(iterations)      
    %%
    fprintf('iteration no. %d/%d \n',k,iterations);
    c=0;
    % Implement first optimal control move and update plant states.
    [x_k1, u] = getTransitionsTrue(x,c,x_w,r0,Kstable);
    [xl_k1, ul] = getTransitions(xl,c,Kstable);
    % shift the output so that it's from the working point perspective
    % setpoint being [0;0;0;0]
    %%%%%%%%%%%%%% DATA %%%%%%%%%%%%%
    X=[x(1:2)-x_w(1:2); u-r0];
    Y=(x_k1-x_w)-(A*(x-x_w)+B*(u-r0));
    q=50; % data window of 50 datapoints
    data=update_data(X,Y,q,k,data);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xo_k1,uo]=getTransitionsLearn(xo,c,Kstable,data);

    his = [x-x_w; u-r0]; % c: decision var; u-r0: delta u;
    hisO=[xo;uo];
    hisL=[xl;ul];
    % Save plant states for display.
    sysHistory = [sysHistory his]; %#ok<*AGROW>
    sysHistoryL = [sysHistoryL hisL]; %#ok<*AGROW>
    sysHistoryO = [sysHistoryO hisO]; %#ok<*AGROW>
    
    
    xo=xo_k1;
    x=x_k1;
    xl=xl_k1;
end
toc

%% PLOT
figure;
subplot(n+m,1,1);
plot(0:iterations,sysHistory(1,:),'Linewidth',1.5); hold on;
plot(0:iterations,sysHistoryL(1,:),'Linewidth',1.5,'LineStyle','--'); hold on;
plot(0:iterations,sysHistoryO(1,:),'Linewidth',1.5,'LineStyle','-.','Color','g');
grid on
xlabel('iterations');
ylabel('x1');
title('mass flow');
subplot(n+m,1,2);
plot(0:iterations,sysHistory(2,:),'Linewidth',1.5); hold on;
plot(0:iterations,sysHistoryL(2,:),'Linewidth',1.5,'LineStyle','--'); hold on;
plot(0:iterations,sysHistoryO(2,:),'Linewidth',1.5,'LineStyle','-.','Color','g');
grid on
xlabel('iterations');
ylabel('x2');
title('pressure rise');
subplot(n+m,1,3);
plot(0:iterations,sysHistory(3,:),'Linewidth',1.5); hold on;
plot(0:iterations,sysHistoryL(3,:),'Linewidth',1.5,'LineStyle','--'); hold on;
plot(0:iterations,sysHistoryO(3,:),'Linewidth',1.5,'LineStyle','-.','Color','g');
grid on
xlabel('iterations');
ylabel('x3');
title('throttle');
subplot(n+m,1,4);
plot(0:iterations,sysHistory(4,:),'Linewidth',1.5); hold on;
plot(0:iterations,sysHistoryL(4,:),'Linewidth',1.5,'LineStyle','--'); hold on;
plot(0:iterations,sysHistoryO(4,:),'Linewidth',1.5,'LineStyle','-.','Color','g');
grid on
xlabel('iterations');
ylabel('x4');
title('throttle rate');
subplot(n+m,1,5);
plot(0:iterations,sysHistory(5,:),'Linewidth',1.5); hold on;
plot(0:iterations,sysHistoryL(5,:),'Linewidth',1.5,'LineStyle','--'); hold on;
plot(0:iterations,sysHistoryO(5,:),'Linewidth',1.5,'LineStyle','-.','Color','g');

grid on
xlabel('iterations');
ylabel('u');
title('Sys input');
legend({'True system', 'Linear system', 'Learned true system'});


figure;
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1.5,'LineStyle','-'); hold on;
plot(sysHistoryL(1,:),sysHistoryL(2,:),'Linewidth',1.5,'LineStyle','--');  hold on;
plot(sysHistoryO(1,:),sysHistoryO(2,:),'Linewidth',1.5,'LineStyle','-.','Color','g');
grid on
xlabel('x1');
ylabel('x2');
title('State space');


function data=update_data(X,Y,q,iter,data_old)
    % Moving window of q datapoints
    if iter<q
    % UPDATE DATA for estimation
        data.X=[data_old.X, X];
        data.Y=[data_old.Y,Y];
    else
        data.X=[data_old.X(:,2:end), X];
        data.Y=[data_old.Y(:,2:end),Y];
    end
end