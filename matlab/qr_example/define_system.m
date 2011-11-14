function [A,B,C,d_0] = define_system(add_altitude_dynamics)
%%
%==========================================================================
% Define A, B, C of discrete time model:
% x+ = Ax + Bu
% y = Cx
%
%==========================================================================

% A, B based on Pat's ME290J project
A_block = [...
    1 0.025 0.0031 0;
    0 1     0.2453 0;
    0 0     0.7969 0.0225;
    0 0    -1.7976 0.9767];

B_block = [0; 0; 0.01; 0.9921];

if add_altitude_dynamics
    m = 1.3; % kg
    dT = 0.025; % s
    
    % add altitude dynamics:
    [Aalt, Balt, Dalt] = altitude_dynamics(m, dT);
    A = blkdiag(A_block, A_block, Aalt);
    B = [B_block zeros(size(B_block)) zeros(size(B_block)); ...
         zeros(size(B_block)) B_block zeros(size(B_block)); ...
         zeros(size(Balt)) zeros(size(Balt)) Balt];
    C = eye(10);
    
    %% Initial disturbance estimate
    d_0 = [0 0 0 0  0 0 0 0 Dalt']';

else
    A = blkdiag(A_block, A_block);
    B = [B_block zeros(size(B_block)); ...
        zeros(size(B_block)) B_block];
    C = eye(8);
    %% Initial disturbance estimate
    d_0 = [0 0 0 0  0 0 0 0]';

end

