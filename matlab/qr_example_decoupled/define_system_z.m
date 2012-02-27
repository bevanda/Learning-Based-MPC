m = 1.3; % kg
dT = 0.025; % s

% Altitude dynamics
g = 9.81; % N/kg
eff_factor = 0.91;

Ac = [0 1; 0 0];
Bc = [0; -eff_factor/m];

dc = [0; g];

%dT = 0.025;

Ad = expm(Ac*dT);

intexpA = [dT dT^2/2; 0 dT];

Bd = intexpA*Bc;

dd = intexpA*dc;

A = Ad;
B = Bd;
C = eye(2);
d_0 = dd;
% Count number of states n, number of inputs m, number of outputs o:
n = size(A,1);
m = size(B,2);
o = size(C,1);


%% Get design parameters:
% Define parameters of MPC procedure

Q = diag([...
    20 ... % z pos
    1  % z vel
    ]); % Cost on states
%Q = Q(1:p,1:p);
R = diag([...
    0.1 ... % thrust cmd
    ]); % Cost on control
R = R(1:m,1:m);

dlqr_controlweight = 8; % cost multiplier for 'baseline' feedback law ** maybe try 5
maxadm_controlweight = 100; % cost multiplier for generating LQR K that is
% used in terminal set computation
max_z = -0.13; % [m] NED convention
min_z = -3.0; % [m] NED convention
max_vz = 5; % [m/s]
max_thrust_cmd = 18; % N
min_thrust_cmd = 0; %N

% Uncertainty bound on x+ (use same for x and y axes)
state_uncert = [...
    0.01 ... % alt [m] % this should be conservative
    0.1 ... % alt rate [m/s] % based on max_thrust*beta_max(7) = 18*0.0035 = 0.0630, and some conservatism
    ]';


%%
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================

temp0 = [max_thrust_cmd];
temp01 = [-min_thrust_cmd];
temp0 = temp0(1:m);
temp01 = temp01(1:m);
F_u = [eye(m); -eye(m)]; h_u = [temp0;temp01];
temp1 = [max_z; max_vz];
temp2 = [-min_z; max_vz];
F_x = [eye(n); -eye(n)]; h_x = [temp1;temp2];
F_g = [eye(n); -eye(n)]; h_g = [state_uncert; state_uncert]; % uncertainty polytope