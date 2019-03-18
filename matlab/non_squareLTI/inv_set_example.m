%% Calculating of the maximal admissible invariant set O_{inf}(0)

A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
% Count number of states n, number of inputs m, number of outputs o:
n = size(A,1);
m = size(B,2);
o = size(C,1);

Ntheta = [1, 0];
% MN = [Mtheta; Ntheta];
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
% Mtheta = [1, 0, 0, 0; 0, 1, 1, -2]';
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);

Q = diag([1,1]);
R = diag([1,1]);
K = -dlqr(A, B, Q, R);
sys1 = LTISystem('A', A+B*K);
umax = [0.3;0.3]; umin = [-0.3;-0.3];
xmax = [5; 5]; xmin = [-5; -5];
R = Polyhedron([K; -K; eye(2); -eye(2)], [umax; -umin; xmax; -xmin]);
Q = Polyhedron.unitBox(2)*4;
P = R&Q; % intersection of P and Q
iset1 = sys1.invariantSet('X', P);
iset.plot();

%% Calculating maximal admissible invariant set for tracking

ALPHA = 0.99; 
L = [-K eye(m)]*[LAMBDA' PSI']';
A_w = [A+B*K, B*L;...
       zeros(n,m), eye(m)];
   
sys2 = LTISystem('A', A_w);
umax = [0.3;0.3]; umin = [-0.3;-0.3];
xmax = [5; 5]; xmin = [-5; -5];
aug_state_constraints =
R = Polyhedron([[eye(2), zeros(2,2)]; [zeros(2,2),-eye(2)]; []], [xmax; -xmin; umax; -umin; ALPHA*]);
Q = Polyhedron.unitBox(4)*4;
P = R&Q; % intersection of P and Q
iset2 = sys2.invariantSet('X', P);
iset2.plot() 