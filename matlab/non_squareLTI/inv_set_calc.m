%% Simple nonsquare discrete LTI system 
% Count number of states n, number of inputs m, number of outputs o:
A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
n = size(A,1);
m = size(B,2);
o = size(C,1);

Q = diag([1,1]);
R = diag([1,1]);
K = -dlqr(A, B, Q, R);
% %% Calculating of the maximal admissible invariant set O_{inf}(0)
% 
% Ntheta = [1, 0];
% % MN = [Mtheta; Ntheta];
% M = [A - eye(n), B, zeros(n,o); ...
%         C, zeros(o,m), -eye(o)];
% Mtheta = null(M);
% % Mtheta = [1, 0, 0, 0; 0, 1, 1, -2]';
% LAMBDA = Mtheta(1:n,:);
% PSI = Mtheta(n+1:n+m,:);
% 
% sys1 = LTISystem('A', A+B*K);
% umax = [0.3;0.3]; umin = [-0.3;-0.3];
% xmax = [5; 5]; xmin = [-5; -5];
% R = Polyhedron([K; -K; eye(2); -eye(2)], [umax; -umin; xmax; -xmin]);
% Q = Polyhedron.unitBox(2)*4;
% P = R&Q; % intersection of P and Q
% iset1 = sys1.invariantSet('X', P);
% iset1.plot();
%% Tracking solution with artificial setpoint

% MN = [Mtheta; 1, 0];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);
%%%%%%%%%%%%%%%%%%%%%
% UNDER DISTURBANCE %
d_0 = [0,0]';
% Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + V*theta, theta in R^m
V_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = V_0(1:n);
PSI_0 = V_0(n+1:n+m);
%%%%%%%%%%%%%%%%%%%%%

L = (PSI - K*LAMBDA);
L0 = (PSI_0 - K*LAMBDA_0); % when being under inital disturbance

umax = [0.3;0.3]; umin = [-0.3;-0.3];
xmax = [5; 5]; xmin = [-5; -5];

F_u = [eye(m); -eye(m)]; h_u = [umax; -umin];
F_x = [eye(n); -eye(n)]; h_x = [xmax; -xmin];

% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);


F = [F_x     zeros(length_Fx, m);...
     F_x*K     F_u*L; ...
     zeros(length_Fx, n) F_x*LAMBDA; ...
     zeros(length_Fu, n) F_u*PSI];

% contract the h values for the artificial steady state by a scalar λ ∈ (0, 1)
lambda=0.99;
h = [...
    h_x; ...
    h_u; ...
    lambda*h_x; ...
    lambda*h_u];
Ak=A+B*K;
Aw = [A+B*K     B*L; ...
      zeros(m,n) ones(m)];
% sys = LTISystem('Ak', Ak);
P = Polyhedron([F_x;F_u], [h_x;h_u]);

% iset = sys.invariantSet('X', P,'maxIterations',10);
% assert(~iset.isEmptySet(),"The invariant set is empty");
iset = P;
for n = 1:50
    iset = iset + Ak^n*P;
end
figure
iset.plot()
P = Polyhedron(F, h);
syss = LTISystem('A', Aw);
iset1 = syss.invariantSet('X', P,'maxIterations',30);
iset1=iset1.projection(1:2);
figure
iset1.plot()
%% INTERSETING calculation for set X_N
% sysStruct.A=A+B*K;
% sysStruct.B=B;
% sysStruct.C=C;
% sysStruct.D=D;
% sysStruct.xmax = xmax;
% sysStruct.xmin = xmin;
% sysStruct.umax = umax;
% sysStruct.umin= umin;
% % sysStruct.Q = Q;
% % sysStruct.R = R;
% system = LTISystem(sysStruct);
% InvSet = system.invariantSet();
% termCons_F=InvSet.A;
% termCons_h=InvSet.b;
% InvSet.plot();