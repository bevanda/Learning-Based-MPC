%==========================================================================
% Define all design parameters
%==========================================================================
function [...
    Q, R, ...
    dlqr_controlweight, maxadm_controlweight, ...
    x_max, ...
    u_max, ...
    state_uncert, ...
    enable_learning, ALPHA, MU_factor, ...
    uncertainty_block] ...
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
   
    dlqr_controlweight = 8; % cost multiplier for 'baseline' feedback law ** maybe try 5
    maxadm_controlweight = 100; % cost multiplier for generating LQR K that is
    % used in terminal set computation
    x_max = 5;
    u_max = 0.3;
    
%% 
% Uncertainty bound on x+ (use same for x and y axes)

uncert = [...
    0.1 ... % mass flow
    0.1 ... % pressure rise
    ];
state_uncert = uncert';


ALPHA = 0.99; % forgetting factor (closer to 1, 'forget' slower)
enable_learning = 1;
MU_factor = 0.001; % weight on prior model ** maybe try 0.01

uncertain_struct_option = 7;

% no longer used:
% uncertainty_upper_bounds = 0.1*ones(p+m+1,1);
% uncertainty_lower_bounds = -0.1*ones(p+m+1,1);
%% WHAT IS THIS UNCERTAINTY STRUCTURE?????????????????????????????????????????????????????????????????????
%??????????????????????????????????????????????/?
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
    
elseif uncertain_struct_option == 13

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
    
    elseif uncertain_struct_option == 7
    % new
    uncertainty_block = [...
        % A       B  d
        0 0 0 0   0  1
        0 0 0 0   0  1
        0 0 0 0   0  1
        0 0 1 1   0  1
        ];
end

