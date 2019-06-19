% Solves the minimization problem with T=10
%  minimize       1/2*integral{t=0 until t=T}(x1^2 + x2^2 + u^2) dt
%  subject to     dot(x1) = (1-x2^2)*x1 - x2 + u,   x1(0)=0, x1(T)=0
%                 dot(x2) = x1,                     x2(0)=1, x2(T)=0
%                 x1(t) >= -0.25, 0<=t<=T
% Joel Andersson, UW Madison 2017
%
% States and control
import casadi.*

x1 = casadi.SX.sym('x1');
x2 = casadi.SX.sym('x2');
u = casadi.SX.sym('u');

% Model equations
x1_dot = (1-x2^2)*x1 - x2 + u;
x2_dot = x1;

% Least squares objective terms
lsq =[x1; x2; u];

% Define the problem structure
ocp = struct('x',[x1; x2],'u',u, 'ode', [x1_dot; x2_dot], 'lsq',lsq);

% Specify problem data
data = struct('T', 10,...
              'x0', [0;1],...
              'xN', [ 0; 0],...
              'x_min', [-0.25; -inf],...
	      'x_max', [ inf;  inf],...
	      'x_guess',  [0; 0],...
              'u_min', -1,...
              'u_max', 1,...
              'u_guess', 0);

% Specify solver options
opts = struct('N', 20,...
	      'verbose', true);

% Create an OCP solver instance
s = GNMS(ocp, data, opts);

% Initializing figure
figure();
clf;
hold on;
grid on;

% Plot solution
x1_plot = plot(s.t, s.sol.x(1,:), 'r--');
x2_plot = plot(s.t, s.sol.x(2,:), 'b-');
u_plot = stairs(s.t, [s.sol.u(1,:) nan], 'g-.');

xlabel('t')
legend('x1','x2','u')
pause(2);

iter=0;
while s.sol.norm_dw>1e-8
     % SQP iteration
     s.sqpstep();
     iter = iter + 1;

     % Update plots
     set(x1_plot, 'Ydata', s.sol.x(1, :));
     set(x2_plot, 'Ydata', s.sol.x(2, :));
     set(u_plot, 'Ydata', s.sol.u(1,:));
     title(sprintf('Iteration %d, |dw| = %g', iter, s.sol.norm_dw))
     pause(2);
end