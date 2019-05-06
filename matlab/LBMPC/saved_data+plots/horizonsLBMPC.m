close all;
addpath('./data/')
%% Load data
LBMPC50r=load('LBMPC_N50_sys_full.mat');
LBMPC60r=load('LBMPC_N40_sys_term.mat');
LBMPC80r=load('LBMPC_N40_sys_full.mat');

%% Plot
iterations=size(LBMPC80r.sysH,2)-1;
Ts=0.01;
t=Ts*(0:iterations);

figure;
subplot(5,1,1);
plot(t, LBMPC50r.sysH(1,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t, LBMPC60r.sysH(1,:),'Linewidth',1.5,'Color','m'); hold on;
plot(t,LBMPC80r.sysH(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(t,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_1');
title('mass flow');

subplot(5,1,2);
plot(t, LBMPC50r.sysH(2,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t, LBMPC60r.sysH(2,:),'Linewidth',1.5,'Color','m'); hold on;
plot(t, LBMPC80r.sysH(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(t, 0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_2');
title('pressure rise');
legend({'LBMPC N=50' 'LBMPC N=60','LBMPC N=80'},'Location','northwest')

subplot(5,1,3);
plot(t, LBMPC50r.sysH(3,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t, LBMPC60r.sysH(3,:),'Linewidth',1.5,'Color','m'); hold on;
plot(t,LBMPC80r.sysH(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(t,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(t, LBMPC50r.sysH(4,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t, LBMPC60r.sysH(4,:),'Linewidth',1.5,'Color','m'); hold on;
plot(t, LBMPC80r.sysH(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(t, 0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(t, LBMPC50r.sysH(5,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t, LBMPC60r.sysH(5,:),'Linewidth',1.5,'Color','m'); hold on;
plot(t,LBMPC80r.sysH(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(t, 0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta u');
title('Sys input');



