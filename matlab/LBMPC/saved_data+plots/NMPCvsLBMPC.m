close all;
addpath('./data/')
%% Load data
LBMPC40=load('LBMPC_N40_sys_full.mat');
LBMPC50=load('LBMPC_N50_sys_full.mat');

NMPC40=load('NMPC_N40_sys_full.mat');
NMPC50=load('NMPC_N50_sys_full.mat');

%% Plot
iterations=size(NMPC50.sysH,2)-1;
Ts=0.01;
t=Ts*(0:iterations);

figure;
subplot(5,1,1);
plot(t, LBMPC40.sysH(1,:),'Linewidth',1.7,'Color','r'); hold on;
plot(t, LBMPC50.sysH(1,:),'Linewidth',1.7,'Color','b'); hold on;
plot(t, NMPC40.sysH(1,:),'Linewidth',1.7,'Color','r','LineStyle',':'); hold on;
plot(t, NMPC50.sysH(1,:),'Linewidth',1.7,'Color','b','LineStyle',':'); hold on;
plot(t,0*(0:iterations),'Linewidth',1.5,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_1');
title('mass flow');


subplot(5,1,2);
plot(t, LBMPC40.sysH(1,:),'Linewidth',1.7,'Color','r'); hold on;
plot(t, LBMPC50.sysH(1,:),'Linewidth',1.7,'Color','b'); hold on;
plot(t, NMPC40.sysH(1,:),'Linewidth',1.7,'Color','r','LineStyle',':'); hold on;
plot(t, NMPC50.sysH(1,:),'Linewidth',1.7,'Color','b','LineStyle',':'); hold on;
plot(t,0*(0:iterations),'Linewidth',1.5,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_2');
title('pressure rise');
legend({'LBMPC N=40', 'LBMPC N=50', 'NMPC N=40','NMPC N=50'},'Location','northwest')

subplot(5,1,3);
plot(t, LBMPC40.sysH(3,:),'Linewidth',1.7,'Color','r'); hold on;
plot(t, LBMPC50.sysH(3,:),'Linewidth',1.7,'Color','b'); hold on;
plot(t, NMPC40.sysH(3,:),'Linewidth',1.7,'Color','r','LineStyle',':'); hold on;
plot(t, NMPC50.sysH(3,:),'Linewidth',1.7,'Color','b','LineStyle',':'); hold on;
plot(t,0*(0:iterations),'Linewidth',1.5,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(t, LBMPC40.sysH(4,:),'Linewidth',1.7,'Color','r'); hold on;
plot(t, LBMPC50.sysH(4,:),'Linewidth',1.7,'Color','b'); hold on;
plot(t, NMPC40.sysH(4,:),'Linewidth',1.7,'Color','r','LineStyle',':'); hold on;
plot(t, NMPC50.sysH(4,:),'Linewidth',1.7,'Color','b','LineStyle',':'); hold on;
plot(t,0*(0:iterations),'Linewidth',1.5,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(t, LBMPC40.sysH(5,:),'Linewidth',1.7,'Color','r'); hold on;
plot(t, LBMPC50.sysH(5,:),'Linewidth',1.7,'Color','b'); hold on;
plot(t, NMPC40.sysH(5,:),'Linewidth',1.7,'Color','r','LineStyle',':'); hold on;
plot(t, NMPC50.sysH(5 ,:),'Linewidth',1.7,'Color','b','LineStyle',':'); hold on;
plot(t,0*(0:iterations),'Linewidth',1.5,'LineStyle',':','Color','k'); 
grid on;
xlabel('time');
ylabel('\delta u');
title('Sys input');

