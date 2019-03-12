clear all; clc; close all;
%% INIT CONTROLLER DESIGN
syms u ... % control input
    x1 ... % mass flow
    x2 ... % pressure rise
    x3 ... % throttle opeining
    x4 ... % throttle opening rate
    xk1 ...
    xk2 ...
    xk3 ...
    xk4 ...
    ;
wn=sqrt(1000); % resonant frequency
zeta=1/sqrt(2); % damping coefficient
beta=1; % constant >0
x2_c=0; % pressure constant
%% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;
%% Continous time state-space model of the Moore-Greitzer compressor model

f1 = x2+x2_c+1+3*(x1/2)-(x1^3/2); % mass flow rate
f2 = (x1+1-x3*sqrt(x2))/(beta^2); % pressure rise rate
f3 = x4; % throttle opening rate
f4 = -wn^2*x3-2*zeta*wn*x4+wn^2*u; % throttle opening acceleration

%% Linearisation around the equilibrium [0.5 1.6875 1.1547 0]'

A = jacobian([f1,f2, f3, f4], [x1, x2, x3, x4]);
B = jacobian([f1,f2, f3, f4], [u]);
% equilibrium params
x1 = 0.5;
x2 = 1.6875;
x3 = 1.1547;
x4 = 0;
init_cond = [x1-0.35, x2-0.4, x3, 0];  % init condition
% print the matrices in the cmd line
A = eval(A);
B = eval(B);
% C = [A(1, :); A(2,:)] % choose f1 and f2 as outputs 
C = eye(4);
D = zeros(4,1);
eig(A)
Kc=-place(A,B,[-22.5,-21.5,-1.0,-0.1])
% Visualise the poles and zeros of the continuous system
[b,a]=ss2tf(A,B,C,D);
sys = tf([b(1,:)],[a]);
% sys2 = tf([b(2,:)],[a]);
% figure;
% pzmap(sys);
% grid on;
% pzmap(sys2);
%% Euler discretisation

Ts = 0.01; % sampling time
xk=[xk1,xk2,xk3,xk4]';
Ad = eye(4)+A*Ts;
Bd = Ts*B;
Cd = eye(4);
Dd = D;
e = eig(Ad);
% figure;
% sys = idss(Ad,Bd,Cd,Dd,'Ts',0.01);
% pzmap(sys);

%% System stabilisation /w feedback matrix K to place poles at ~ p=[0.75, 0.78, 0.98, 0.99]

p=[0.99, 0.98, 0.78, 0.75]; % desired poles of the open-loop system, while still being stable close to existing ones
[K,prec,message] = place(Ad,Bd,p); %nominal feedback matrix
% K=[-3.0741 2.0957 0.1197 -0.0090]; %nominal feedback matrix from the LBMPC paper
AK = Ad-Bd*K;
e = eig(AK);
xkk = (AK)*xk;
% figure;
% sys = idss(AK,zeros(4,1),Cd,Dd,'Ts',0.01);
% pzmap(sys);
%% Init for data in simulink
X =[0 0 0];
Y=[0 0 0 0];