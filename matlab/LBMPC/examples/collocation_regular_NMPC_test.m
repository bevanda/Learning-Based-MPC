%
% regular NMPC collocation for Moore-Greitzer comperssor model
clearvars; clc; 
addpath('../utilities');
addpath('../models/'); 
addpath('../functions/'); 

import casadi.*

% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0 collocation_points(d, 'legendre')];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);
end

% Parameter initialization
disp('Initializing parameters...');
% Generate the DLTI verison of the continouous Moore-Greitzer Compressor 
% model (MGCM)
[A,B,~,~,Ts]=mgcmDLTI();
% system dimensions
n = size(A,1); % num of states
m = size(B,2); % num of inputs
o = size(C,1); % num of outputs
% Obtaining all needed matrices for the optimal control problem (OCP)
[Kstabil,Klqr,Q,R,P,T,Mtheta,LAMBDA,PSI,LAMBDA_0,PSI_0]=matOCP(A,B,C,n,m,o);


% Time horizon
Nt = 0.5;

% eqilibilium point
x_eq = [0.500000000000000;1.68750000000000;1.15470000000000;0];
u_eq = 1.15470000000000;
x_init = [0.150000000000000;1.28750000000000;1.15470000000000;0];

% Constraints of the compressor model
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;

umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];

% Declare model variables
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x = [x1; x2; x3; x4];
u = SX.sym('u');

% Model equations
xdot = [-x2+1+3*(x1/2)-(x1^3/2);...       % mass flow rate 
        (x1+1-x3*sqrt(x2));...            % pressure rise rate 
        x4;...                                % throttle opening rate
        -1000*x3-2*sqrt(500)*x4+1000*u];    % throttle opening accelerat
L = (x1-x_eq(1))^2 + (x2-x_eq(2))^2 + (x3-x_eq(3))^2 + (x4-x_eq(4))^2  + (u-u_eq)^2;  

% Continuous time dynamics (TRUE)
f = Function('f', {x, u}, {xdot, L});

% Control discretization
N = 50; % number of control intervals
h = Nt/N;

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
Xk = MX.sym('X0', n);
w = {w{:}, Xk};
lbw = [lbw; x_init];
ubw = [ubw; x_init];
w0 = [w0; x_init];

% Formulate the NLP
for k=0:N-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw; umin];
    ubw = [ubw;  umax];
    w0 = [w0;  x_init(3)];

    % State at collocation points
    Xkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], n);
        w = {w{:}, Xkj{j}};
        lbw = [lbw; xmin];
        ubw = [ubw;  xmax];
        w0 = [w0; zeros(n,1)];
    end

    % Loop over collocation points
    Xk_end = D(1)*Xk;
    for j=1:d
       % Expression for the state derivative at the collocation point
       xp = C(1,j+1)*Xk;
       for r=1:d
           xp = xp + C(r+1,j+1)*Xkj{r};
       end

       % Append collocation equations
%        [fj, qj] = f(Xkj{j},Uk); % true
       [fj, qj] = lin_f(Xkj{j},Uk); % nominal
       g = {g{:}, h*fj - xp};
       lbg = [lbg; zeros(n,1)];
       ubg = [ubg; zeros(n,1)];

       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkj{j};

       % Add contribution to quadrature function
       J = J + B(j+1)*qj*h;
    end

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], n);
    w = {w{:}, Xk};
    lbw = [lbw; xmin];
    ubw = [ubw;  xmax];
    w0 = [w0; zeros(n,1)];

    % Add equality constraint
    g = {g{:}, Xk_end-Xk};
    lbg = [lbg; zeros(n,1)];
    ubg = [ubg; zeros(n,1)];
end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

% Plot the solution
x1_opt = w_opt(1:(m+n)+n*d:end);
x2_opt = w_opt(2:(m+n)+n*d:end);
x3_opt = w_opt(3:(m+n)+n*d:end);
x4_opt = w_opt(4:(m+n)+n*d:end);
u_opt = w_opt(5:(m+n)+n*d:end);
tgrid = linspace(0, Nt, N+1);
plotRESPONSE([x1_opt';x2_opt';x3_opt';x4_opt'; [u_eq; u_opt]'],tgrid,n,m);
