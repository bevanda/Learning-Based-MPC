% clear workspace, close open figures
clearvars
close all
clc

%% import CasADi
import casadi.*
addpath('/home/petar/acados/external/casadi-matlabR2014b-v3.4.0/');

%%

% A non-square DLTI system - sampled double integrator
A = [1 1; 0 1];
B = [0.0 0.5; 1.0 0.5];
C = [1 0];
% system dimensions
n = size(A,1); % state dimension
m = size(B,2); % input dimension 
o = size(C,1);
Q = eye(n);
R = eye(m);

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


%%

% Number of MPC iterations
mpciterations = 100;

% Horizon (discrete)
N = 3;

% initial conditions
t_0 = 0.0;
% The initial conditions
x_init = [0;...
         -2];
%setpoint
xs = [4.95;...
      0.0];
%% Define a nominal feedback policy K and corresponding terminal cost
K = -dlqr(A, B, Q, R); % 'baseline' stabilizing feedback gain
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
% terminal steady state cost
T = 100*P;
% stage cost
Q = eye(n);
R = eye(m);

% Terminal set and cost
K_loc = K;

% Initial guess for input
u0 = zeros(m*N,1); % start inputs
theta0 = zeros(m,1); % start param values

%initial state constraint: use LB, UB
%input constraints (actually put in nl_cons())
lb=-inf*ones(m*(N+1),1);
ub=+inf*ones(m*(N+1),1);
%nonlinear constraints (both inequality and equality constraints)
con_lb=zeros(m*(N+1),1); % to formulate as in fmincon cieq <= 0
con_ub=inf*ones(m*(N+1),1); % to formulate as in fmincon cieq <= 0
%make symbolic
y=MX.sym('y',N*m+m);
global X; % system state
X = x_init;
obj=cost_fun(N, y, X, xs, A,B,Q, R, P, T, m, LAMBDA, PSI);
con=nl_cons(N, y, X, m, A, B, F_run,h_run, F_term, h_term);
nlp = struct('x', c, 'f', obj, 'g', con);
solver = nlpsol('solver', 'ipopt', nlp); %,'file_print_level',5

% Set variables for output
t = [];
x = [];
u = [];
theta =[];
% ellipsoids toolbox needed (Matlab central)
%E = ellipsoid(x_eq, alpha*inv(P));

f1 = figure(1); hold on
set(f1,'PaperPositionMode','auto')
set(f1,'Units','pixels')
% plot terminal set
%plot(E,'r'), axis equal, grid on


% Print Header
fprintf('   k  |      u(k)        theta(k)    Time \n');
fprintf('---------------------------------------------------\n');

% initilization of measured values
tmeasure = t_0;
xmeasure = x_init;

%% simulation
for ii = 1:mpciterations % maximal number of iterations
    
    
    % Set initial guess and initial constraint
    fprintf("Iterations : %d\n",k);
    xs = set_ref(k);
    y_init=[u0;theta0];
    
    t_Start = tic;
    % calling 'solver' get the solution
    res = solver('x0' , y_init,... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg', con_lb,...           % lower bound on g
             'ubg', con_ub);             % upper bound on g
    y_OL=full(res.x); 
    u_OL=y_OL(1:m*(N));
    theta_OL=y_OL(m*N+1:end);
    t_Elapsed = toc( t_Start );    
    %%    
 
    % Store closed loop data
    t = [ t, tmeasure ];
    u_new=u_OL(1:m);
    u = [ u, u_new ];
    theta_new=u_OL(1:m);
    theta = [ theta, theta_new ];
    xmeasure = dynamics(xmeasure,u_new,A,B);
    x = [ x, xmeasure];
    X = xmeasure;
    % Update closed-loop system (apply first control move to system)
    tmeasure = tmeasure + delta;
    % Compute initial guess for next time step, based on terminal LQR controller (K_loc)
    u0 = [u_OL(m+1:end); K_loc*x_measure];
    theta0 = theta_OL;
    %%
    % Print numbers
    fprintf(' %3d  |  %+11.6f %+11.6f %+11.6f %+11.6f  %+6.3f\n', ii, u(1,end),u(2,end),...
            x(1,end), x(2,end),t_Elapsed);
    
    %plot predicted and closed-loop state trajetories    
    f1 = figure(1);
    plot(x(1,:),x(2,:),'b'), grid on, hold on,
    plot(x_OL(1:n:n*(N+1)),x_OL(n:n:n*(N+1)),'g')
    plot(x(1,:),x(2,:),'ob')
    xlabel('x(1)')
    ylabel('x(2)')
    drawnow
  
end
% plot
figure(2)
stairs(t,u(1)); hold on;
stairs(t,u(2)); grid on;
figure(3)
plot(t,x(1,:),'b');
grid on;
hold on;
plot(t,theta(1,:),'b');


