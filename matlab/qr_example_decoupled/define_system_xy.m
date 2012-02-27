% System definition
A = [...
    1 0.025 0.0031 0;
    0 1     0.2453 0;
    0 0     0.7969 0.0225;
    0 0    -1.7976 0.9767];

B = [0; 0; 0.01; 0.9921];
C = eye(4);
d_0 = [0 0 0 0]'; % Initial disturbance estimate
    
% Count number of states n, number of inputs m, number of outputs o:
n = size(A,1);
m = size(B,2);
o = size(C,1);


%% Get design parameters:
% Define parameters of MPC procedure

Q = diag([...
    20 ...  % x pos
    1 ... % x vel
    0.01 ... % pitch
    0.01 ... % pitch rate
    ]); % Cost on states
%Q = Q(1:p,1:p);
R = diag([...
    1 ... % roll cmd
    ]); % Cost on control
R = R(1:m,1:m);

dlqr_controlweight = 8; % cost multiplier for 'baseline' feedback law ** maybe try 5
maxadm_controlweight = 100; % cost multiplier for generating LQR K that is
% used in terminal set computation
max_x = 2; % [m]
max_vx = 10; % [m/s]
max_pitch_cmd = deg2rad(35); % [rad]
max_pitch = pi/2; %pi/4; % [rad]
max_pitch_rate = 4*pi; % [rad/s]

% Uncertainty bound on x+ (use same for x and y axes)
state_uncert = [...
    0.01 ... % pos [m]
    0.1 ... % vel [m/s]
    0.01 ... % att [rad]
    0.1 % att rate [rad/s]
    ]';


%%
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x.  Also define model uncertainty as a F_g*x <= h_g
%==========================================================================

temp0 = [max_pitch_cmd];
temp01 = [max_pitch_cmd];
temp0 = temp0(1:m);
temp01 = temp01(1:m);
F_u = [eye(m); -eye(m)]; h_u = [temp0;temp01];
temp1 = [max_x; max_vx; max_pitch; max_pitch_rate];
temp2 = [max_x; max_vx; max_pitch; max_pitch_rate];
F_x = [eye(n); -eye(n)]; h_x = [temp1;temp2];
F_g = [eye(n); -eye(n)]; h_g = [state_uncert; state_uncert]; % uncertainty polytope