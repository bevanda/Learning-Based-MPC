%% Calculating of the maximal admissible invariant set O_{inf}(0)
A = [1, 1; 0, 1];
B = [0.0, 0.5; 1.0, 0.5];
C = [1 0];
D = [0 0];
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
[P, e, K] = dare(A,B,Q,R);
sys1 = LTISystem('A', A-B*K);
umax = [0.3;0.3]; umin = [-0.3;-0.3];
xmax = [5; 5]; xmin = [-5; -5];
% R = Polyhedron([K; -K; eye(2); -eye(2)], [umax; -umin; xmax; -xmin]);
% Q = Polyhedron.unitBox(2)*4;
% P = R&Q; % intersection of P and Q
% iset1 = sys1.invariantSet('X', P);
% iset1.plot();

%% Calculating maximal admissible invariant set for tracking

% ALPHA = 0.99; 
% L = [-K eye(m)]*[LAMBDA' PSI']';
% A_w = [A+B*K, B*L;...
%        zeros(n,m), eye(m)];
%    
% sys2 = LTISystem('A', A_w);
% umax = [0.3;0.3]; umin = [-0.3;-0.3];
% xmax = [5; 5]; xmin = [-5; -5];
% R = Polyhedron([[eye(2), zeros(2,2)]; [zeros(2,2),-eye(2)]; []], [xmax; -xmin; umax; -umin]);
% Q = Polyhedron.unitBox(4)*4;
% P = R&Q; % intersection of P and Q
% iset2 = sys2.invariantSet('X', P);
% iset2.plot() 

%% Calculating maximal admissible control invariant set
sysStruct.A=A-B*K;
sysStruct.B=B;
sysStruct.C=C;
sysStruct.D=D;
sysStruct.xmax = xmax;
sysStruct.xmin = xmin;
sysStruct.umax = umax;
sysStruct.umin= umin;
% sysStruct.Q = Q;
% sysStruct.R = R;
system = LTISystem(sysStruct);
InvSet = system.invariantSet();
termCons_F=InvSet.A;
termCons_h=InvSet.b;
InvSet.plot();

%% Calculating MACI set with parameter dependence 

ALPHA = 0.99; 

Mtheta =[LAMBDA' PSI']';

L = [K eye(m)]*Mtheta;

sysStruct.A=[A-B*K ,     B*L;...
            zeros(n,m), eye(m)]^3;
sysStruct.B=zeros(4,2);
sysStruct.C=zeros(1,4);
sysStruct.D=zeros(1,2);

sysStruct.xmax = [xmax;inf;inf]*ALPHA;
sysStruct.xmin = [xmin;-inf;-inf]*ALPHA;
sysStruct.umax = umax*ALPHA;
sysStruct.umin= umin*ALPHA;
% sysStruct.umax = umax*ALPHA;
% sysStruct.umin= umin*ALPHA;
% sysStruct.x.penalty.weight=Q;
% sysStruct.u.penalty.weight=R;

system = LTISystem(sysStruct);
InvSet2 = system.invariantSet(); % InvSet2 is a polyhaeder
% extracting H-representation
term_F=InvSet2.A;
term_h=InvSet2.b;
% InvSet2.plot();

% project the 4D case to a 2D one
MAI=projection(InvSet2,1:2); % Maximal Admissible Invariant set
plot(MAI);