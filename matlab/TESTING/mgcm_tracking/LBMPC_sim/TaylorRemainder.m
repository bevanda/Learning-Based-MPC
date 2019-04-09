close all; clear vars;
%% INIT CONTROLLER DESIGN
syms u ... % control input
    x1 ... % mass flow
    x2 ... % pressure rise
    x3 ... % throttle opeining
    x4 ... % throttle opening rate
    ;
wn=sqrt(1000); % resonant frequency
zeta=1/sqrt(2); % damping coefficient
beta=1; % constant >0https://www.google.com/search?client=ubuntu&channel=fs&q=yonotope&ie=utf-8&oe=utf-8
x2_c=0; % pressure constant
%% Constraints
mflow_min=0; mflow_max=1;
prise_min=1.1875; prise_max=2.1875;
throttle_min=0.1547; throttle_max=2.1547;
throttle_rate_min=-20; throttle_rate_max=20;
u_min=0.1547;u_max=2.1547;
%% Continous time state-space model of the Moore-Greitzer compressor model

f1 = -x2+x2_c+1+3*(x1/2)-(x1^3/2); % mass flow rate
f2 = (x1+1-x3*sqrt(x2))/(beta^2); % pressure rise rate
f3 = x4; % throttle opening rate
f4 = -wn^2*x3-2*zeta*wn*x4+wn^2*u; % throttle opening acceleration
ff= [f1;f2;f3;f4];
xx=[x1;x2;x3;x4];
%% Linearisation around the equilibrium [0.5 1.6875 1.1547 0]'
A = jacobian([f1,f2, f3, f4], [x1, x2, x3, x4]);
B = jacobian([f1,f2, f3, f4], [u]);

% equilibrium params / 1-order  taylor expansion point   
x10 = 0.5;
x20 = 1.6875;
x30 = 1.1547;
x40 = 0;
u=x30;
equili = [x10; x20; x30; x40];
init_cond = [x10-0.35, x20-0.4, x30, 0];  % init condition

%% Linearization vs true sys plot for f1
figure;
% 1st order  taylor expansion point  | in matlab is order 2
T1=taylor(f1,xx(1:2),'Order',2,'ExpansionPoint', equili(1:2));
ezl1=ezplot(T1,[mflow_min,mflow_max,prise_min,prise_max]);
hold on;
T1=simplify(T1);
% constraints as plot limits 
eyt1=ezplot(f1,[mflow_min,mflow_max,prise_min,prise_max]);
legend('Linarization','True')
set(ezl1,'color',[1 0 0])

%% Linearization vs ture sys plot for f2
figure;
% 1st order  taylor expansion point  | in matlab is order 2
f2=(x1+1-x30*sqrt(x2))/(beta^2); % to set r to working point value to have 2 var for plotting
T2=taylor(f2,xx(1:2),'Order',2,'ExpansionPoint', equili(1:2));
ezl2=ezplot(T2,[mflow_min,mflow_max,prise_min,prise_max]);
hold on;
T2=simplify(T2);
% constraints as plot limits 
ezt2=ezplot(f2,[mflow_min,mflow_max,prise_min,prise_max]);
legend('Linarization','True')
set(ezl2,'color',[1 0 0])

%% Calculating the Taylor approx error
% X1
E1=(f1-T1);
E1=simplify(E1^2);
figure;
fplot(E1);
% constrained optimization /w symbolic variables
objfun1=matlabFunction(-E1,'vars',x1);
[x1, fval1]=fmincon(objfun1,[0],[],[],[],[],[mflow_min],[mflow_max],[],[]);
-fval1
% X2
E2=(f2-T2);
E2=simplify(E2^2);
figure;
fplot(E2);
% constrained optimization /w symbolic variables
objfun2=matlabFunction(-E2,'vars',x2);
[x2, fval2]=fmincon(objfun2,[0],[],[],[],[],[prise_min],[prise_max],[],[]);
-fval2