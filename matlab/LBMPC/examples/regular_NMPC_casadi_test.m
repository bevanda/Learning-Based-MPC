% clear workspace, close open figures
clearvars
close all
clc

%% import CasADi
import casadi.*
addpath('/home/pbevanda/casadi-linux-matlabR2014b-v3.4.5/');

%%

% system dimensions
n = 4; % state dimension
m = 1; % input dimension 

% eqilibilium point
x_eq = [0.500000000000000;1.68750000000000;1.15470000000000;0];
u_eq = 1.15470000000000;

mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;

% Number of MPC iterations
mpciterations = 500;

% Horizon (continuous)
T = 0.5;

% sampling time (Discretization steps)
delta = 0.01;

% Horizon (discrete)
N = T/delta;

% initial conditions
t_0 = 0.0;
x_init = [0.150000000000000;1.28750000000000;1.15470000000000;0];

% stage cost
Q = eye(n);
R = eye(m);

% Terminal set and cost
K_loc = [-3.07418713694076,2.09578024408799,0.119436236659434,-0.00894688869207627];
P = [2379.75328799635,-1408.52636933972,586.834405598336,2.54452796979169;
    -1408.52636933972,1020.35695741999,-486.636840230204,-2.14515528702556;
    586.834405598336,-486.636840230204,376.289369319808,1.71242577067112;
    2.54452796979169,-2.14515528702556,1.71242577067112,1.01279182904087];
%alpha = 0.0250; % Caution! Computation may last several minutes per iteration!
alpha = 0.7000; % obtained via the alternative way

% Initial guess for input
u0 = u_eq*ones(m*N,1);
% Initial gues for states by simulation
x0=zeros(n*(N+1),1);
x0(1:n) = x_init;
for k=1:N
     x0(n*k+1:n*(k+1)) = dynamic(delta,x0(n*(k-1)+1:n*k), u0(k));
end

%initial state constraint: use LB, UB
%input constraints
lb=[-inf*ones(n*(N+1),1);(u_eq-1)*ones(m*N,1)];
ub=[+inf*ones(n*(N+1),1);(u_eq+1)*ones(m*N,1)];
lb(1:n)=x_init;
ub(1:n)=x_init;
%nonlinear constraints (both inequality and equality constraints)

xmax = [zeros(n,1); mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [zeros(n,1); mflow_min; prise_min; throttle_min; throttle_rate_min];

con_lb=repmat(xmin,N,1);
con_ub=repmat(xmax,N,1);

%make symbolic
y=MX.sym('y',N*m+(N+1)*n);
obj=costfunction(N, y, x_eq, u_eq, Q, R, P,n,m,delta);
con=nonlinearconstraints(N, delta, y, x_eq, P, alpha,n,m);
nlp = struct('x', y, 'f', obj, 'g', con);
solver = nlpsol('solver', 'ipopt', nlp); %,'file_print_level',5

% Set variables for output
t = [];
x = [];
u = [];

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

% simulation
for ii = 1:mpciterations % maximal number of iterations
    
    
    % Set initial guess and initial constraint
    beq=xmeasure;
    y_init=[x0;u0];
    
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
    u_OL=y_OL(n*(N+1)+1:end);
    t_Elapsed = toc( t_Start );    
    %%    
 
    % Store closed loop data
    t = [ t, tmeasure ];
    x = [ x, xmeasure ];
    u = [ u, u_OL(1:m) ];
    
    % Update closed-loop system (apply first control move to system)
    xmeasure = x_OL(n+1:2*n);
    tmeasure = tmeasure + delta;
        
    % Compute initial guess for next time step, based on terminal LQR controller (K_loc)
    u0 = [u_OL(m+1:end); K_loc*x_OL(end-n+1:end)];
    x0 = [x_OL(n+1:end); dynamic(delta, x_OL(end-n+1:end), u0(end-m+1:end))];
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
%%
figure(2)
stairs(t,u);
figure(3)
stairs(x(1,:)); hold on;
stairs(x(2,:));
%%
function xdot = system(x, u)
    % Systemn dynamics
    xdot = [ -x(2)+1+3*(x(1)/2)-(x(1)^3/2);... % mass flow rate 
            (x(1)+1-x(3)*sqrt(x(2)));... % pressure rise rate 
            x(4);... % throttle opening rate
            -1000*x(3)-2*sqrt(500)*x(4)+1000*u]; % throttle opening acceleration];
    
end

function cost = costfunction(N, y, x_eq, u_eq, Q, R, P,n,m,delta)
    % Formulate the cost function to be minimized
    
    cost = 0;
    x=y(1:n*(N+1));
    u=y(n*(N+1)+1:end);
    
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N
        x_k=x(n*(k-1)+1:n*k);
        u_k=u(m*(k-1)+1:m*k);
        cost = cost + delta*runningcosts(x_k, u_k, x_eq, u_eq, Q, R);
    end
    cost = cost + terminalcosts( x(n*N+1:n*(N+1)), x_eq, P);
    
end

function cost = runningcosts(x, u, x_eq, u_eq, Q, R)
    % Provide the running cost   
    cost = (x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
    
end


  function [con] = nonlinearconstraints(N, delta, y, x_eq, P, alpha,n,m) 
   % Introduce the nonlinear constraints also for the terminal state
   
   x=y(1:n*(N+1));
   u=y(n*(N+1)+1:end);
   con = [];
   %con_ub = [];
   %con_lb = [];
   % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        u_k=u((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_new - dynamic(delta, x_k, u_k);
        ceqnew2 = x_new;
        con = [con; ceqnew; ceqnew2];
        
        %nonlinear constraints on state and input could be included here
    end
   %
   %terminal constraint
%    [cnew] = terminalconstraints( x(n*N+1:n*(N+1)), x_eq, P, alpha);
%     con = [con; cnew];    
end

function cost = terminalcosts(x, x_eq, P)
    % Introduce the terminal cost
    cost = (x-x_eq)'*P*(x-x_eq);
end


function [con] = terminalconstraints(x, x_eq, P, alpha)
    % Introduce the terminal constraint
    con   = (x-x_eq)'*P*(x-x_eq) - alpha;
end

function x_new=dynamic(delta,x,u)
    %use Ruku4 for discretization
    k1=system(x,u);
    k2=system(x+delta/2*k1,u);
    k3=system(x+delta/2*k2,u);
    k4=system(x+delta*k3,u);
    x_new=x+delta/6*(k1+2*k2+2*k3+k4);
end
