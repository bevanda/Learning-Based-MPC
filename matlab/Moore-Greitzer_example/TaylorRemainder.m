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

A = jacobian([f1,f2, f3, f4], [x1, x2, x3, x4])
B = jacobian([f1,f2, f3, f4], [u])


% equilibrium params / 1-order  taylor expansion point   
x1 = 0.5;
x2 = 1.6875;
x3 = 1.1547;
x4 = 0;
u=x3;
equili = [x1; x2; x3; x4];
init_cond = [x1-0.35, x2-0.4, x3, 0];  % init condition
T1=taylor(ff,xx,'Order',1)
T1=simplify(T1);
%%
syms x
g = exp(x*sin(x));
t = taylor(g, 'ExpansionPoint', 1.5, 'Order', 2);
t = simplify(t);
size(char(t))
xd = 1:0.05:3;
yd = subs(g,x,xd);
fplot(t, [1, 3])
hold on
plot(xd, yd, 'r-.')
title('Taylor approximation vs. actual function')
legend('Taylor','Function')
