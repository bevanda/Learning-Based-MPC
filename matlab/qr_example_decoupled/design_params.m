%==========================================================================
% Define all design parameters
%==========================================================================
function [...
    N, Q, R, ...
    dlqr_controlweight, maxadm_controlweight, ...
    max_x, max_vx, ...
    max_pitch_cmd, ...
    max_pitch, max_pitch_rate,  ...
    state_uncert] = design_params(p,m)

% Define parameters of MPC procedure
N = 15; % MPC horizon

weights_select = 0; % 0 = ball catching, 1 = testing, 2 = learning demo (back & forth)

    Q = diag([...
        20 ...  % x pos
        1 ... % x vel
        0.01 ... % pitch
        0.01 ... % pitch rate
        ]); % Cost on states
    Q = Q(1:p,1:p);
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


%%


% Uncertainty bound on x+ (use same for x and y axes)
state_uncert = [...
    0.01 ... % pos [m]
    0.1 ... % vel [m/s]
    0.01 ... % att [rad]
    0.1 % att rate [rad/s]
    ]';
