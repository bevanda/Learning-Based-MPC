%==========================================================================
% Define all design parameters
%==========================================================================
function [...
    Q, R, ...
    x_max, x_min, ...
    u_max, u_min, ...
    state_uncert, ...
    ALPHA, MU_factor] ...
    = design_params(p,m) % p -state dim, m -input dim

% Define parameters of MPC procedure
%N = 15; % MPC horizon

    Q = diag([...
        1 ...  % mass flow
        1 ... % pressure rise
        ]); % Cost on states
    Q = Q(1:p,1:p);
    R = diag([...
        1 ...
        1 ...
        ]); % Cost on control
    R = R(1:m,1:m);
    % used in terminal set computation
    x_max = 5; x_min = -5;
    u_min = 0.3; u_max = -0.3;
    
%% 
% Uncertainty bound on x+ (use same for x and y axes)

uncert = [...
    0.0 ... % x1
    0.0 ... % x2
    ];
state_uncert = uncert';


ALPHA = 0.99; % forgetting factor (closer to 1, 'forget' slower)

MU_factor = 0.001; % weight on prior model ** maybe try 0.01



end
