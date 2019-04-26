function [Ad,Bd,Cd,Dd,Ts]=mgcmDLTI()
%--------------------------------------------------------------------------
% Discrete time linear time invariant model of the Moore-Greitzer
% Compressor Model
%--------------------------------------------------------------------------
syms u ... % control input
    x1 ... % mass flow
    x2 ... % pressure rise
    x3 ... % throttle opeining
    x4 ... % throttle opening rate
    ;
wn=sqrt(1000); % resonant frequency
zeta=1/sqrt(2); % damping coefficient
beta=1; % constant >0
x2_c=0; % pressure constant

% Continous time state-space model of the Moore-Greitzer compressor model
f1 = -x2+x2_c+1+3*(x1/2)-(x1^3/2); % mass flow rate 
f2 = (x1+1-x3*sqrt(x2))/(beta^2); % pressure rise rate 
f3 = x4; % throttle opening rate
f4 = -wn^2*x3-2*zeta*wn*x4+wn^2*u; % throttle opening acceleration

% Linearisation around the equilibrium [0.5 1.6875 1.1547 0]'
A = jacobian([f1,f2, f3, f4], [x1, x2, x3, x4]);
B = jacobian([f1,f2, f3, f4], [u]);
% equilibrium values
x1 = 0.5;
x2 = 1.6875;
x3 = 1.1547;
x4 = 0;
A = eval(A);
B = eval(B);
C = eye(4);
D = zeros(4,1);
n=size(A,2);
% Exact discretisation
Ts = 0.01; % sampling time
Ad = expm(A*Ts);
Bd = (Ad-eye(n))*inv(A)*B;
Cd=C;
Dd=D;

end
