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

% Uncertainty bound obtained form the maximal linearization error
% Lagrange Error Bound + estimation error toleranace
state_uncert = [0.02;5e-04;0;0];


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
 F_x_d,h_x_d]... % uncertainty ineq
    =getCONSPOLY(...
    xmax,xmin,umax,umin,state_uncert,...
    x_wp,u_wp,m,n,...
    A,B,Q,R,LAMBDA,PSI,LAMBDA_0,PSI_0);

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

% Initial guess for input
u0 = u_eq*ones(m*N,1);
theta0=zeros(m,1);
% Initial gues for states by simulation
x0=zeros(n*(N+1),1);
x0(1:n) = x_init;
for k=1:N
     x0(n*k+1:n*(k+1)) = x_eq + nominal_dynamics(x0(n*(k-1)+1:n*k)-x_eq, u0(k)-u_eq,A,B);
end
xl0 = x0; % learned model = nominal at the beginning

%input constraints
lb=[-inf*ones(2*n*(N+1),1);-inf*ones(N*m+m,1)];
ub=[+inf*ones(2*n*(N+1),1);+inf*ones(N*m+m,1)];
lb(1:n)=x_init;
ub(1:n)=x_init;
lb((N+1)*n+1:(N+1)*n+n)=x_init;
ub((N+1)*n+1:(N+1)*n+n)=x_init;
%nonlinear constraints (both inequality and equality constraints)

xmin = [zeros(n,1); zeros(n,1); -inf*ones(numel(h_x)+numel(h_u),1)];
xmax = [zeros(n,1); zeros(n,1); zeros(numel(h_x)+numel(h_u),1)]; 

con_lb=[-inf*ones(numel(h_w_N) + numel(h_x_d),1); repmat(xmin,N,1)];
con_ub=[zeros(numel(h_w_N) + numel(h_x_d),1); repmat(xmax,N,1)];

y=SX.sym('y',2*(N+1)*n+N*m+m);
q=100; % moving window of q datapoints 
d = SX.sym('d',8,q); %data as parameter
obj=costfunction(N, y, x_eq, u_eq,  Q, R, P,T, LAMBDA, PSI, n, m, delta);
con=nonlinearconstraints(N,A, B, d,x_eq,u_eq,y,n,m,...
    F_x,h_x, F_u,h_u, F_w_N, h_w_N, F_x_d, h_x_d);

nlp = struct('x', y, 'f', obj, 'g', con,'p',d);
solver = nlpsol('solver', 'ipopt', nlp); %,'file_print_level',5

% ellipsoids toolbox needed (Matlab central)
%E = ellipsoid(x_eq, alpha*inv(P));

%% simulation

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
xlmeasure = x_init;

% Set logging vectors
t = [];
x = [];
xl = [];
u = [];
theta = [];
art_ref = [];
solve_times=[];
data = zeros(8,q); 
data(8,1)=1; % take into account init data

for iter = 1:mpciterations % maximal number of iterations

    % Set initial guess and initial constraint
    y_init=[xl0;x0;u0;theta0];
    
    t_Start = tic;
    lb(1:n)=xmeasure;
    ub(1:n)=xmeasure;
    lb(n*(N+1)+1:n*(N+1)+n)=xmeasure;
    ub(n*(N+1)+1:n*(N+1)+n)=xmeasure;
    
    res = solver('x0' , y_init,...   % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg', con_lb,...       % lower bound on g
             'ubg', con_ub,...       % upper bound on g
             'p', data);             % parameters
    y_OL=full(res.x); 
    xl_OL=y_OL(1:n*(N+1));
    x_OL=y_OL(n*(N+1)+1:2*n*(N+1));
    u_OL=y_OL(2*n*(N+1)+1:end-m);
    theta_OL=y_OL(end-m+1:end);
    art_ref_OL=[LAMBDA*y_OL(end-m+1:end); PSI*y_OL(end-m+1:end)];
    t_Elapsed = toc( t_Start );   
    u_in = u_OL(1:m); % get first input
    
    %%    
    solve_times = [solve_times,t_Elapsed];
    % Store closed loop data
    t = [ t, tmeasure ];
    x = [ x, xmeasure ];
    xl = [ xl, xlmeasure ];
    u = [ u, u_in ];
    theta = [theta, theta_OL];
    art_ref = [art_ref, art_ref_OL];
    
    % Update closed-loop system (apply first control move to system)
    xmeasure1 = dynamic(delta,xmeasure,u_OL(1:m)); % real state measurement
    xlmeasure = x_eq + learned_dynamics(xmeasure-x_eq, u_in - u_eq, A, B, data);
    
    u_o = u_OL(1:m);
     % data acquisition 
    X=[xmeasure(1:2)-x_eq(1:2); u_o-u_eq]; %[δx1;δx2;δu]
    Y=xmeasure1-(x_eq+A*(xmeasure-x_eq)+ B*(u_o-u_eq)); %[δx_true-δx_nominal]
    data = get_data(X,Y,q,iter,data); % update data      
    
    % Compute initial guess for next time step, based on terminal LQR controller (Kstabil)
    u0 = [u_OL(m+1:end); u_eq + Kstabil*(x_OL(end-n+1:end))];
    x0 = [x_OL(n+1:end); x_eq + nominal_dynamics(x_OL(end-n+1:end)-x_eq, u_OL(end-m+1:end)- u_eq, A, B)]; %% learned ynamics
    xl0 = [xl_OL(n+1:end); x_eq + learned_dynamics(xl_OL(end-n+1:end)-x_eq, u_OL(end-m+1:end) - u_eq, A, B, data)]; %% learned dynamics
    theta0 = theta_OL;
    art_ref_OL=[LAMBDA*theta_OL; PSI*theta_OL];

    tmeasure = tmeasure + delta;
    xmeasure = xmeasure1; % UPDATE STATE MEASURE
   
    
    %%
    % Print numbers
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f  %+6.3f\n', iter, u(end),...
            x(1,end), x(2,end),t_Elapsed);
    %plot predicted and closed-loop state trajetories    
    f1 = figure(1);
    plot(x(1,:),x(2,:),'b'), grid on, hold on,
    plot(x_OL(1:n:n*(N+1)),x_OL(2:n:n*(N+1)),'g')
    plot(xl_OL(1:n:n*(N+1)),xl_OL(2:n:n*(N+1)),'r')
    plot(x(1,:),x(2,:),'ob')
    xlabel('x_{1}')
    ylabel('x_{2}')
    drawnow
  
end

%% plotting
plot_RESPONSE([x;u], art_ref+[x_eq;u_eq], t, n, m);
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

function cost = costfunction(N, y, x_eq, u_eq, Q, R, P,T,LAMBDA,PSI, n,m,delta)
    % Formulate the cost function to be minimized
    cost = 0;
    x=y(1:n*(N+1));
    u=y(2*n*(N+1)+1:end-m);
    theta=y(end-m+1:end);
    x_art=LAMBDA*theta;
    u_art=PSI*theta;
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N
        x_k=x(n*(k-1)+1:n*k);
        u_k=u(m*(k-1)+1:m*k);
        cost = cost + delta*runningcosts(x_k - x_eq, u_k - u_eq, x_art, u_art, Q, R);
    end
    cost = cost + terminalcosts(x(n*N+1:n*(N+1)), x_eq, x_art + x_eq, P,T);
    
end
    
function cost = runningcosts(x, u, x_art, u_art, Q, R)
    % Provide the running cost   
    cost = (x-x_art)'*Q*(x-x_art) + (u-u_art)'*R*(u-u_art);
    
end

function cost = terminalcosts(x,x_eq, x_art, P,T)
    % Introduce the terminal cost
    cost = (x-x_art)'*P*(x-x_art)+(x_eq-x_art)'*T*(x_eq-x_art);
end

function [con] = nonlinearconstraints(N, A,B,data,x_eq,u_eq, y,n,m,...
      state_F,state_h, in_F,in_h, F_w_N, h_w_N, F_x_d, h_x_d) 
   % Introduce the nonlinear constraints also for the terminal state
   
   xl=y(1:n*(N+1));
   x=y(n*(N+1)+1:2*n*(N+1));
   u=y(2*n*(N+1)+1:end-m);
   theta=y(end-m+1:end);
   con = [];
   %con_ub = [];
   %con_lb = [];
   % constraints along prediction horizon
    for k=1:N
        xl_k=xl((k-1)*n+1:k*n);
        xl_new=xl(k*n+1:(k+1)*n);   
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n); 
        u_k=u((k-1)*m+1:k*m);
        if k == 1
            cieq_run = F_x_d*(x_new-x_eq)-h_x_d;
            cieq_T = F_w_N*[(x_new-x_eq);theta]-h_w_N;
            con = [con; cieq_run; cieq_T]; %#ok<*AGROW>
        end
        % dynamic constraint
        ceqnew1=xl_new - (x_eq + learned_dynamics(xl_k-x_eq, u_k-u_eq,A,B,data));
        ceqnew2=x_new - (x_eq + nominal_dynamics(x_k-x_eq, u_k-u_eq,A,B));
        con = [con; ceqnew1; ceqnew2];
        % other constraints
        cieq_run1 = state_F*(x_new-x_eq)-state_h;
        cieq_run2 = in_F*(u_k-u_eq)-in_h;
        con  = [con; cieq_run1; cieq_run2]; %#ok<*AGROW>
        % nonlinear constraints on state and input could be included here
    end
  end

 function xk1=nominal_dynamics(xk, uk, A, B)
    xk1 = A*xk + B*uk;
 end
 
 function [xk, u] = learned_dynamics(x, u, A,B, data)
%% Discrete-time learned dynamic model 
% Inputs:
%   x: states at time k
%   u: input at time k
%   data: struct with (X, Y) obervations for estimation
%
% Outputs:
%   xk: states at time k+1
%
%% Discrete time state-space model of the non-square LTI system for tracking
xk = A*x + B*u + casadiL2NW(x,u,data);
 end
 


function x_new=dynamic(delta,x,u)
    %use Ruku4 for discretization
    k1=system(x,u);
    k2=system(x+delta/2*k1,u);
    k3=system(x+delta/2*k2,u);
    k4=system(x+delta*k3,u);
    x_new=x+delta/6*(k1+2*k2+2*k3+k4);
end
