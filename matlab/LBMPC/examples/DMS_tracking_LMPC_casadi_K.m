% clear workspace, close open figures
clearvars;
close all;
clc;
%% import CasADi
import casadi.*

addpath('../utilities');
addpath('../models/'); 
addpath('../functions/'); 
%% Parameter initialization
disp('Initializing parameters...');
% Generate the DLTI verison of the continouous Moore-Greitzer Compressor 
% model (MGCM)
[A,B,C,D,Ts]=mgcmDLTI();
% system dimensions
n = size(A,1); % num of states
m = size(B,2); % num of inputs
o = size(C,1); % num of outputs

% Obtaining all needed matrices for the optimal control problem (OCP)
[Kstabil,Klqr,Q,R,P,T,Mtheta,LAMBDA,PSI,LAMBDA_0,PSI_0]=matOCP(A,B,C,n,m,o);

%% Problem preparation

disp('Setting up the OPC...');

% Constraints of the compressor model
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;

umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];

% Uncertainty bound 
state_uncert = [0;0;0;0]; % no uncertainty

% The initial conditions w.r.t. to the linearisation/working point
x_wp_init = [-0.35;...
            -0.4;...
            0.0;...
            0.0];
% setpoint w.r.t. to the linearisation/working point
x_wp_ref = [0.0;...
            0.0;...
            0.0;...
            0.0];
 
% Working point (wp)
x_wp = [0.5;...
        1.6875;...
        1.1547;...
        0.0];
u_wp = x_wp(3);

[F_x,h_x, ... % nominal state ineq constraints 
 F_u,h_u,...  % nominal input ineq constraints 
 F_w_N,h_w_N,... % terminal extended state ineq constraints 
]... 
    =getCONS(...
    xmax,xmin,umax,umin,...
    x_wp,u_wp,m,n,...
    A,B,Kstabil,LAMBDA,PSI,LAMBDA_0,PSI_0);

%% Setting up the OCP


% eqilibilium point
x_eq = [0.500000000000000;1.68750000000000;1.15470000000000;0];
u_eq = 1.15470000000000;


% Number of MPC iterations
mpciterations = 500;

% Time horizon (continuous)
N_t = 1.0;

% sampling time (Discretization steps)
delta = 0.01;

% Horizon (discrete)
N = N_t/delta;

% initial conditions
t_0 = 0.0;
x_init = [0.150000000000000;1.28750000000000;1.15470000000000;0];

% Terminal set and cost
K_loc = Kstabil;

% Initial guess for input
c0 = zeros(N*m,1);
u0 = u_eq*ones(m*N,1);
theta0=zeros(m,1);
% Initial gues for states by simulation
x0=zeros(n*(N+1),1);
x0(1:n) = x_init;
for k=1:N
     x0(n*k+1:n*(k+1)) = x_eq+...
         nominal_dynamics(x0(n*(k-1)+1:n*k)-x_eq, u0(k)-u_eq,A,B); 
end

%initial state constraint: use LB, UB
%input constraints
lb=[-inf*ones(n*(N+1),1);-inf*ones(2*N*m+m,1)];
ub=[+inf*ones(n*(N+1),1);+inf*ones(2*N*m+m,1)];

lb(1:n)=x_init;
ub(1:n)=x_init;
% lb(n*(N+1)+1:n*(N+1)+m)=zeros(m,1);
% ub(n*(N+1)+1:n*(N+1)+m)=zeros(m,1);

%nonlinear constraints (both inequality and equality constraints)
run_cons_ub = [zeros(m,1); zeros(n,1); zeros(numel(h_x)+numel(h_u),1)]; %[u_e, x_e, term_cons]
run_cons_lb = [zeros(m,1); zeros(n,1); -inf*ones(numel(h_x)+numel(h_u),1)]; %[u_e, x_e, term_cons]

con_lb=[repmat(run_cons_lb,N,1);-inf*ones(numel(h_w_N),1)];
con_ub=[repmat(run_cons_ub,N,1);zeros(numel(h_w_N),1)];

%make symbolic
y=SX.sym('y',(N+1)*n + 2*N*m + m);
obj=costfunction(N, y, Kstabil, x_eq, u_eq,  Q, R, P,T, LAMBDA,PSI, n,m,delta);
con=nonlinearconstraints(N,A,B, Kstabil, delta,x_eq,u_eq,y,n,m,F_x,h_x, F_u,h_u, F_w_N, h_w_N);

nlp = struct('x', y, 'f', obj, 'g', con);
solver = nlpsol('solver', 'ipopt', nlp); %,'file_print_level',5

% Set variables for output
t = [];
x = [];
c = [];
u = [];
theta = [];
art_ref =[];
solve_times=[];
% ellipsoids toolbox needed (Matlab central)
%E = ellipsoid(x_eq, alpha*inv(P));

f1 = figure(1); hold on
set(f1,'PaperPositionMode','auto')
set(f1,'Units','pixels')
% plot terminal set
%plot(E,'r'), axis equal, grid on

% Print Header
fprintf('   k  |      u(k)        x(1)        x(2)     Time \n');
fprintf('---------------------------------------------------\n');

% initilization of measured values
tmeasure = t_0;
xmeasure = x_init;
%% simulation
for ii = 1:mpciterations % maximal number of iterations
    
    
    % Set initial guess and initial constraint
    beq=xmeasure;
    y_init=[x0;c0;u0;theta0];
    
    t_Start = tic;
    lb(1:n)=xmeasure;
    ub(1:n)=xmeasure;
    res = solver('x0' , y_init,... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg', con_lb,...           % lower bound on g
             'ubg', con_ub);             % upper bound on g
    y_OL=full(res.x); 
    x_OL=y_OL(1:n*(N+1));
    c_OL=y_OL(n*(N+1)+1:n*(N+1)+N*m);
    u_OL=y_OL(n*(N+1)+N*m+1:end-m);
    theta_OL=y_OL(end-m+1:end);
    art_ref_OL=[LAMBDA*y_OL(end-m+1:end);PSI*y_OL(end-m+1:end)];
    t_Elapsed = toc( t_Start );   
    solve_times = [solve_times,t_Elapsed];
    %%    
    % Store closed loop data u_OL(1:m)
    t = [ t, tmeasure ];
    x = [ x, xmeasure ];
    c = [ c, c_OL(1:m)];
    u = [ u, u_OL(1:m) ];
    theta = [theta, theta_OL];
    art_ref = [art_ref, art_ref_OL];
    
    % Update closed-loop system (apply first control move to system)
    xmeasure = dynamic(delta, xmeasure, u_OL(1:m)); %simulate 
    tmeasure = tmeasure + delta;
    % Compute initial guess for next time step, based on terminal LQR controller (K_loc)
    c0 = zeros(N*m,1);
    u0 = [u_OL(m+1:end); K_loc*x_OL(end-n-m+1:end-m)];
    x0 = [x_OL(n+1:end); (x_eq+nominal_dynamics(x_OL(end-n-m+1:end-m)-x_eq, u0(end-m-m+1:end-m)-u_eq,A,B))];
    theta0 = theta_OL;
    art_ref_OL=[LAMBDA*theta_OL; PSI*theta_OL];
    %%
    % Print numbers
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f  %+6.3f\n', ii, u(end),...
            x(1,end), x(2,end),t_Elapsed);
    
    %plot predicted and closed-loop state trajetories    
    f1 = figure(1);
    plot(x(1,:),x(2,:),'b'), grid on, hold on,
    plot(x_OL(1:n:n*(N+1)),x_OL(2:n:n*(N+1)),'g')
    plot(x(1,:),x(2,:),'ob')
    xlabel('x(1)')
    ylabel('x(2)')
    drawnow
  
end

%% plotting
xl=x;
plot_RESPONSE([xl;u], art_ref+[x_eq;u_eq], t, n, m)
%%
fprintf('Total solving time: %6.3fs \n', sum(solve_times));
figure; histfit(solve_times);

%% Help funcs

function xdot = system(x, u)
    % Systemn dynamics
    xdot = [-x(2)+1+3*(x(1)/2)-(x(1)^3/2);...       % mass flow rate 
            (x(1)+1-x(3)*sqrt(x(2)));...            % pressure rise rate 
            x(4);...                                % throttle opening rate
            -1000*x(3)-2*sqrt(500)*x(4)+1000*u];    % throttle opening acceleration
end

function cost = costfunction(N, y, K, x_eq, u_eq, Q, R, P,T,LAMBDA,PSI, n,m,delta)
    % Formulate the cost function to be minimized
    cost = 0;
    x=y(1:n*(N+1));
%     c=y(n*(N+1)+1:n*(N+1)+N*m);
    u=y(n*(N+1)+N*m+1:end-m);
    theta=y(end-m+1:end);
    x_art=LAMBDA*theta;
    u_art=PSI*theta;
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N
        x_k=x(n*(k-1)+1:n*k);
        u_k=u(m*(k-1)+1:m*k);
%         u_k=u((k-1)*m+1:k*m)+K*(x_k-x_eq)+u_eq; % u_hat = K*x+c constraint
        cost = cost + delta*runningcosts(x_k, u_k,x_art+x_eq,u_art+u_eq, Q, R);
    end
    cost = cost + terminalcosts(x(n*N+1:n*(N+1)), x_eq, x_art+x_eq, P,T);
    
end

function cost = runningcosts(x, u, x_art, u_art, Q, R)
    % Provide the running cost   
    cost = (x-x_art)'*Q*(x-x_art) + (u-u_art)'*R*(u-u_art);
    
end

function cost = terminalcosts(x,x_eq, x_art, P,T)
    % Introduce the terminal cost
    cost = (x-x_art)'*P*(x-x_art)+(x_eq-x_art)'*T*(x_eq-x_art);
end


function [con] = nonlinearconstraints(N,A,B, K,delta,x_eq,u_eq, y,n,m,state_F,state_h, in_F,in_h, F_w_N, h_w_N) 
   % Introduce the nonlinear constraints also for the terminal state
   
    x=y(1:n*(N+1));
    c=y(n*(N+1)+1:n*(N+1)+N*m);
    u=y(n*(N+1)+N*m+1:end-m);
    theta=y(end-m+1:end);
    con = [];
    %con_ub = [];
    %con_lb = [];
    % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        u_k=u((k-1)*m+1:k*m);
        c_k=c((k-1)*m+1:k*m);
    %         u_k=u((k-1)*m+1:k*m)+K*(x_k-x_eq)+u_eq; % u_hat = K*x+c constraint
        % dynamic constraint
        % init condition 
        if k == 1
            ceqnew=x_new;
        else
            ceqnew=x_new - (x_eq+nominal_dynamics(x_k-x_eq, u_k-u_eq,A,B));
        end
        ceqnew=x_new - (x_eq+nominal_dynamics(x_k-x_eq, u_k-u_eq,A,B));
        con = [con; ceqnew];
        ceqnew2= u_k - (c_k + u_eq + K*(x_k-x_eq)); % u_hat = K*x+c constraint
        con = [con; ceqnew2];
        % other constraints
        cieq_run1 = state_F*(x_new-x_eq)-state_h;
        cieq_run2 = in_F*(u_k-u_eq)-in_h;
        con  = [con; cieq_run1; cieq_run2]; %#ok<*AGROW>
        % nonlinear constraints on state and input could be included here
    end
    %
    %terminal constraint
    cieq_T = F_w_N*[(x_new-x_eq);theta]-h_w_N; % x-x_eq to be w.r.t. to the eq point
    con = [con; cieq_T];    
    end

function xk1=nominal_dynamics(xk, uk, A, B)
    xk1 = A*xk + B*uk;
 end

function [w] = disturb(w_max,w_min)
    w = rand(4, 1).*(w_max - w_min)+w_min;
end

function x_new=dynamic(delta,x,u)
    %use Ruku4 for discretization
    k1=system(x,u);
    k2=system(x+delta/2*k1,u);
    k3=system(x+delta/2*k2,u);
    k4=system(x+delta*k3,u);
    x_new=x+delta/6*(k1+2*k2+2*k3+k4);
end
