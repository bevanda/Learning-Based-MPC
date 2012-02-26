function [Ad, Bd, dd] = altitude_dynamics(m, dT)
% Altitude dynamics
g = 9.81; % N/kg
%m = 1.3; % kg
eff_factor = 0.91;

Ac = [0 1; 0 0];
Bc = [0; -eff_factor/m];

dc = [0; g];

%dT = 0.025;

Ad = expm(Ac*dT);

intexpA = [dT dT^2/2; 0 dT];

Bd = intexpA*Bc;

dd = intexpA*dc;