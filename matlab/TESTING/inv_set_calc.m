%% Simple nonsquare discrete LTI system 
% Count number of states n, number of inputs m, number of outputs o:
A = [1.01126321746508,-0.0100340214950357,6.46038913508018e-05,1.93716902346107e-07; ...
    0.0100340214950357,0.995515380253533,-0.0127681799951143,-5.57226765949308e-05; ...
    0,0,0.957038195891878,0.00792982548734094; ...
    0,0,-7.92982548734093,0.602405619103784];
B = [4.95338239742896e-07; ...
    -0.000193159646826652; ...
    0.0429618041081219; ...
    7.92982548734093];
C = [1,0,0,0;...
    0,1,0,0;...
    0,0,1,0;...
    0,0,0,1];
% D = [0;0;0;0];
n = size(A,1);
m = size(B,2);
o = size(C,1);

Q = eye(n);
R = eye(m);
%==========================================================================
% Define a nominal feedback policy K and corresponding terminal cost
% 'baseline' stabilizing feedback law
K = -dlqr(A, B, Q, R);
% Terminal cost chosen as solution to DARE
P = dare(A+B*K, B, Q, R);
% terminal steady state cost
T = 1000;
%% Calculating of the maximal admissible invariant set O_{inf}(0)

%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================
% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;

umax = u_max; umin = u_min;
xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];
%  Shift the constraints for the linearised model for the value of the
%  working point
x_w = [0.5;...
    1.6875;...
    1.1547;...
    0.0];
r0 = x_w(3);

F_u = [eye(m); -eye(m)]; h_u = [umax-r0; -umin+r0];
F_x = [eye(n); -eye(n)]; h_x = [xmax-x_w; -xmin+x_w];

Aw = [A+B*K     B*(PSI - K*LAMBDA); ...
      zeros(m,n) ones(m)];
sys = LTISystem('A', Aw);
%==========================================================================
% Compute maximally invariant set
%==========================================================================

disp('Computing and simplifying terminal set...');
F_w = [F_x zeros(length_Fx, m);
    zeros(length_Fx, n) F_x*LAMBDA; ...
    F_u*K, F_u*(PSI - K*LAMBDA); ...
    zeros(length_Fu, n) F_u*PSI];

lambda=0.99; % λ ∈ (0, 1), λ can be chosen arbitrarily close to 1, the obtained
% invariant set can be used as a reliable polyhedral approximation to the maximal invariant set 
h_w = [...
    h_x; ...
    (h_x - F_x*LAMBDA_0)*lambda; ...
    h_u - F_u*(PSI_0 - K*LAMBDA_0); ...
    (h_u - F_u*PSI_0)]*lambda;

P = Polyhedron(F_w, h_w);
iset = sys.invariantSet('X', P);
assert(~iset.isEmptySet(),"The invariant set is empty");
isetX = projection(iset,1:2);
isetX.plot()

%% INTERSETING calculation for set X_N
sysStruct.A=A+B*K;
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