%% TRACKING piecewise constant REFERENCE MPC example
clearvars;
close all;

%% Parameters
% Set the prediction horizon:
N = 3;

% The initial conditions
x = [0;...
    -2];
%setpoint
xs = [4.95;...
      0.0];
options = optimoptions('fmincon','Algorithm','sqp','Display','final');
% Simulation length (iterations)
iterations = 90;
% Optimization variable bounds not constraints per say
u_min =-0.3; u_max= 0.3;
LB = ones(2*N,1)*u_min;
UB = ones(2*N,1)*u_max;

A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
n = size(A,1);
m = size(B,2);
o = size(C,1);

%%%%%%%%%%%%%%%%%%%%%
Q = diag([1,1]);
R = diag([1,1]);

u0 = zeros(m*N,1); % start inputs
opt_var = u0;

%%
%==========================================================================
% Define a nominal feedback policy K and corresponding terminal cost

% 'baseline' stabilizing feedback law
K = -dlqr(A, B, Q, R);
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
%%
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================

umax = [0.3;0.3]; umin = [-0.3;-0.3];
xmax = [5; 5]; xmin = [-5; -5];

F_u = [eye(m); -eye(m)]; h_u = [umax; -umin];
F_x = [eye(n); -eye(n)]; h_x = [xmax; -xmin];

% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);


%%
%==========================================================================
% Compute maximally invariant set
%==========================================================================

disp('Computing and simplifying terminal set...');
F_w = [F_x ; ...
       F_u*K];
% contract the h values for the artificial steady state by a scalar λ ∈ (0, 1)
lambda=0.99;
h_w = [h_x;...
       h_u];

F_w_N0 = F_w; h_w_N0 = h_w;

% Simplify the constraints
term_poly = polytope(F_w_N0, h_w_N0);
[F_w_N, h_w_N] = double(term_poly);
%     term_poly = Polyhedron(F_w_N0, h_w_N0); 
%     F_w_N = term_poly.A; % Inequality description { x | H*[x; -1] <= 0 }   
%     h_w_N = term_poly.b; % Inequality description { x | A*x <= b }

disp('Terminal set Polyhedron:');
term_poly
MAI=projection(term_poly,1:n); % Maximal Admissible Invariant set rpojected on X
plot(MAI);



%% Cost Calculation
% Start simulation
sysHistory = [x;u0(1:2,1)];

true_refHistory = xs;

for k = 1:(iterations)
    xs = set_ref(k);
    COSTFUN = @(var) costFunction(reshape(var,m,N),x,xs,N,reshape(var(1:m),m,1),P,K);
    CONSFUN = @(var) constraintsFunction(reshape(var,m,N),x,N,K,F_w_N,h_w_N);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],LB,UB,CONSFUN,options);    
    c = reshape(opt_var(1:m),m,1);
    % Implement first optimal control move and update plant states.
    x= getTransitions(x, c); 
    his=[x; c];
    % Save plant states for display.
    sysHistory = [sysHistory his]; 
    true_refHistory = [true_refHistory xs];
end


%% Plot

figure;
subplot(3,1,1);
plot(0:iterations,sysHistory(1,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x1');
title('x1');
subplot(3,1,2);
plot(0:iterations,sysHistory(2,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('x2');
title('x2');
subplot(3,1,3);
plot(0:iterations,sysHistory(3,:),0:iterations,sysHistory(4,:),'Linewidth',1);
grid on
xlabel('iterations');
ylabel('u');
title('inputs');

figure;
plot_refs=plot(0:iterations,true_refHistory(1,:),0:iterations,sysHistory(1,:),'Linewidth',1.5);
grid on
xlabel('iterations');
% ylabel('references');
title('true reference vs state response');
legend({'real reference', 'state response'},'Location','northeast')
plot_refs(1).LineStyle='--';
plot_refs(1).Color='Red';
plot_refs(2).Color='Blue';


figure;
plot(sysHistory(1,:),sysHistory(2,:),'Linewidth',1.5,'Marker','.');
grid on
xlabel('x1');
ylabel('x2');
title('State space');

%% Helper functions

%set reference depending on the iteration
function [xs] = set_ref(ct)
    if ct <=30 
        xs=[4.95;0];
    elseif ct > 30 && ct <= 60
        xs=[-5.5;0];
    else
        xs=[2;0];
    end
end