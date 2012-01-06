%==========================================================================
% Define all design parameters
%==========================================================================
function [...
    N, Q, R, ...
    dlqr_controlweight, maxadm_controlweight, ...
    max_x, max_y, max_vx, max_vy, ...
    max_z, min_z, max_vz, ...
    max_pitch_cmd, max_roll_cmd, ...
    max_thrust_cmd, min_thrust_cmd, ...
    max_pitch, max_pitch_rate, max_roll, max_roll_rate, ...
    state_uncert, ...
    enable_learning, ALPHA, MU_factor, ...
    uncertainty_block] = design_params(p,m)

% Define parameters of MPC procedure
N = 60; % MPC horizon

weights_select = 0; % 0 = ball catching, 1 = testing, 2 = learning demo (back & forth)

if weights_select == 0
    Q = diag([...
        20 ...  % x pos
        1 ... % x vel
        0.01 ... % pitch
        0.01 ... % pitch rate
        20 ...  % y pos
        1 ... % y vel
        0.01 ... % roll
        0.01 ... % roll rate
        20 ... % z pos
        1  % z vel
        ]); % Cost on states
    Q = Q(1:p,1:p);
    R = diag([...
        1 ... % roll cmd
        1 ... % pitch cmd
        0.1 ... % thrust cmd
        ]); % Cost on control
    R = R(1:m,1:m);
    
    dlqr_controlweight = 8; % cost multiplier for 'baseline' feedback law ** maybe try 5
    maxadm_controlweight = 100; % cost multiplier for generating LQR K that is
    % used in terminal set computation
    max_x = 2; % [m]
    max_y = 3; % [m]
    max_vx = 10; % [m/s]
    max_vy = 10; % [m/s]
    max_z = -0.13; % [m] NED convention
    min_z = -3.0; % [m] NED convention
    max_vz = 5; % [m/s]
    max_pitch_cmd = deg2rad(35); % [rad]
    max_roll_cmd = deg2rad(35); % [rad]
    max_thrust_cmd = 18; % N
    min_thrust_cmd = 0; %N
    max_pitch = pi/2; %pi/4; % [rad]
    max_pitch_rate = 4*pi; % [rad/s]
    max_roll = pi/2; %pi/4; % [rad]
    max_roll_rate = 4*pi; % [rad/s]


elseif weights_select == 1
    % This is what we tried on Aug 23, 2011:
    
    Q = diag([...
        10 ...  % x pos
        1 ... % x vel
        0.01 ... % pitch
        0.01 ... % pitch rate
        10 ...  % y pos
        1 ... % y vel
        0.01 ... % roll
        0.01 ... % roll rate
        10 ... % z pos
        0.01  % z vel
        ]); % Cost on states
    Q = Q(1:p,1:p);
    R = diag([...
        1 ... % roll cmd
        1 ... % pitch cmd
        0.1 ... % thrust cmd
        ]); % Cost on control
    R = R(1:m,1:m);
    
    dlqr_controlweight = 2; % cost multiplier for 'baseline' feedback law ** maybe try 5
    maxadm_controlweight = 100; % cost multiplier for generating LQR K that is
    % used in terminal set computation
    max_x = 2; % [m]
    max_y = 3; % [m]
    max_vx = 2; % [m/s]
    max_vy = 2; % [m/s]
    max_z = 3; % [m]
    max_vz = 5; % [m/s]
    max_pitch_cmd = deg2rad(20); % [rad]
    max_roll_cmd = deg2rad(20); % [rad]
    max_thrust_cmd = 18; % N
    min_thrust_cmd = 0; %N
    max_pitch = pi/2; %pi/4; % [rad]
    max_pitch_rate = 4*pi; % [rad/s]
    max_roll = pi/2; %pi/4; % [rad]
    max_roll_rate = 4*pi; % [rad/s]
    
elseif weights_select == 2
    % Attempts to de-tune Aug 23, 2011:
    Q = diag([...
        60 ...  % x pos
        1 ... % x vel
        0.01 ... % pitch
        0.01 ... % pitch rate
        60 ...  % y pos
        1 ... % y vel
        0.01 ... % roll
        0.01 ... % roll rate
        60 ... % z pos
        1  % z vel
        ]); % Cost on states
    Q = Q(1:p,1:p);
    R = diag([...
        1 ... % roll cmd
        1 ... % pitch cmd
        0.1 ... % thrust cmd
        ]); % Cost on control
    R = R(1:m,1:m);
    
    dlqr_controlweight = 8; % cost multiplier for 'baseline' feedback law ** maybe try 5
    maxadm_controlweight = 100; % cost multiplier for generating LQR K that is
    % used in terminal set computation
    max_x = 2; % [m]
    max_y = 3; % [m]
    max_vx = 10; % [m/s]
    max_vy = 10; % [m/s]
    max_z = 3; % [m]
    max_vz = 5; % [m/s]
    max_pitch_cmd = deg2rad(35); % [rad]
    max_roll_cmd = deg2rad(35); % [rad]
    max_thrust_cmd = 18; % N
    min_thrust_cmd = 0; %N
    max_pitch = pi/2; %pi/4; % [rad]
    max_pitch_rate = 4*pi; % [rad/s]
    max_roll = pi/2; %pi/4; % [rad]
    max_roll_rate = 4*pi; % [rad/s]

end

%%


% Uncertainty bound on x+ (use same for x and y axes)
axis_uncert = [...
    0.01 ... % pos [m]
    0.1 ... % vel [m/s]
    0.01 ... % att [rad]
    0.1 % att rate [rad/s]
    ];
vert_uncert = [...
    0.01 ... % alt [m] % this should be conservative
    0.1 ... % alt rate [m/s] % based on max_thrust*beta_max(7) = 18*0.0035 = 0.0630, and some conservatism
    ];
if p == 8
    state_uncert = [axis_uncert axis_uncert]';
elseif p == 10
    state_uncert = [axis_uncert axis_uncert vert_uncert]';
else
    assert(0); % shouldn't happen
end

ALPHA = 0.99; % forgetting factor (closer to 1, 'forget' slower)
enable_learning = 1;
MU_factor = 0.001; % weight on prior model ** maybe try 0.01

uncertain_struct_option = 13;

% no longer used:
% uncertainty_upper_bounds = 0.1*ones(p+m+1,1);
% uncertainty_lower_bounds = -0.1*ones(p+m+1,1);

% Uncertainty structure:
uncertainty_block = zeros(p, p + m + 1); % has shape of [A B d]

% By convention, uncertain_struct_option < 10 is for p = 8...
if p == 10
    assert(uncertain_struct_option >= 10);
end

if uncertain_struct_option == 0
    % original
    uncertainty_block = [...
        %      A           B   d
        0 0 0 0  0 0 0 0  0 0  1
        0 0 0 0  0 0 0 0  0 0  1
        0 0 0 0  0 0 0 0  0 0  1
        0 0 1 1  0 0 0 0  1 0  1
        
        0 0 0 0  0 0 0 0  0 0  1
        0 0 0 0  0 0 0 0  0 0  1
        0 0 0 0  0 0 0 0  0 0  1
        0 0 0 0  0 0 1 1  0 1  1
        ];
elseif uncertain_struct_option == 1
    % new
    uncertainty_block = [...
        %      A           B   d
        0 0 0 0  0 0 0 0  0 0  0
        0 0 0 0  0 0 0 0  0 0  1
        0 0 0 0  0 0 0 0  0 0  0
        0 0 1 1  0 0 0 0  1 0  1
        
        0 0 0 0  0 0 0 0  0 0  0
        0 0 0 0  0 0 0 0  0 0  1
        0 0 0 0  0 0 0 0  0 0  0
        0 0 0 0  0 0 1 1  0 1  1
        ];
elseif uncertain_struct_option == 2
    % new
    uncertainty_block = [...
        %      A           B   d
        0 0 0 0  0 0 0 0  0 0  0
        0 0 1 0  0 0 0 0  0 0  1
        0 0 1 1  0 0 0 0  1 0  0
        0 0 1 1  0 0 0 0  1 0  1
        
        0 0 0 0  0 0 0 0  0 0  0
        0 0 0 0  0 0 1 0  0 0  1
        0 0 0 0  0 0 1 1  0 1  0
        0 0 0 0  0 0 1 1  0 1  1
        ];
    % options for altitude dynamics included:
elseif uncertain_struct_option == 11
    % new
    uncertainty_block = [...
        %      A                B     d
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 0 0  0 0  0 0 0  1
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 1 1  0 0 0 0  0 0  1 0 0  1
        
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 0 0  0 0  0 0 0  1
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 1 1  0 0  0 1 0  1
        
        0 0 0 0  0 0 0 0  0 0  0 0 1  1
        0 0 0 0  0 0 0 0  0 0  0 0 1  1
        ];
elseif uncertain_struct_option == 12
    % new
    uncertainty_block = [...
        %      A                B     d
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 0 0  0 0  0 0 0  1
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 1 1  0 0 0 0  0 0  1 0 0  1
        
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 0 0  0 0  0 0 0  1
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 1 1  0 0  0 1 0  1
        
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 0 0  0 0  0 0 1  1
        ];
    
elseif uncertain_struct_option == 13
    % new
    uncertainty_block = [...
        %      A                B     d
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 0 0  0 0  0 0 0  1
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 1 1  0 0 0 0  0 0  1 0 0  1
        
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 0 0  0 0  0 0 0  1
        0 0 0 0  0 0 0 0  0 0  0 0 0  0
        0 0 0 0  0 0 1 1  0 0  0 1 0  1
        
        0 0 0 0  0 0 0 0  0 0  0 0 0  1
        0 0 0 0  0 0 0 0  0 0  0 0 1  1
        ];
end

