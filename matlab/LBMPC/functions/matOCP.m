function [Ks,Klqr,Q,R,P,T,Mtheta,LAMBDA,PSI,LAMBDA_0,PSI_0]=matOCP(A,B,C,n,m,o)
%--------------------------------------------------------------------------
% Obtaining matrices needed for the optimal control problem
%--------------------------------------------------------------------------

% Define a stabilizing nominal feedback policy Ks
p=[0.75, 0.78, 0.98, 0.99]; 
[K,~,~] = place(A,B,p); %nominal feedback matrix
Ks=-K;

% Generate steady-state parameterization
M = [A - eye(n), B, zeros(n,o); ...
        C, zeros(o,m), -eye(o)];
Mtheta = null(M);
LAMBDA = Mtheta(1:n,:);
PSI = Mtheta(n+1:n+m,:);

d_0 = [0,0,0,0]'; % inital disturbance guess
% Solutions of M*[x;u;y] = [-d;0] are of the form M\[-d;0] + M*theta, 
% theta in R^m
M_0 = M\[-d_0; zeros(o,1)];
LAMBDA_0 = M_0(1:n,:);
PSI_0 = M_0(n+1:n+m,:);

% Define a nominal feedback policy K and corresponding terminal cost
% 'baseline' stabilizing feedback law
Q = eye(n); R = eye(m); % state and control penalties
Klqr = -dlqr(A, B, Q, R);

P = dare(A+B*Ks, B, Q, R); % terminal cost chosen as solution to DARE
T = 1000; % terminal steady state cost

end