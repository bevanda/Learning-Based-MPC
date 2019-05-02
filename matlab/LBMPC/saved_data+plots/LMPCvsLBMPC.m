close all;
addpath('./data/')
%% Load data
LMPC50=load('LMPC_N50_sys_full.mat');
LBMPC50=load('LBMPC_N50_sys_full.mat');

%% Plot
iterations=size(LBMPC50.sysH,2)-1;
Ts=0.01;
t=Ts*(0:iterations);

figure;
subplot(5,1,1);
plot(t, LMPC50.sysH(1,:),'Linewidth',1.5,'Color','r'); hold on;
plot(t, LBMPC50.sysH(1,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time ');
ylabel('\delta x_1');
title('mass flow');

subplot(5,1,2);
plot(t, LMPC50.sysH(2,:),'Linewidth',1.5,'Color','r'); hold on;
plot(t, LBMPC50.sysH(2,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time ');
ylabel('\delta x_2');
title('pressure rise');
legend({'LMPC N=50','LBMPC N=50'},'Location','northwest')

subplot(5,1,3);
plot(t, LMPC50.sysH(3,:),'Linewidth',1.5,'Color','r'); hold on;
plot(t, LBMPC50.sysH(3,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time ');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(t, LMPC50.sysH(4,:),'Linewidth',1.5,'Color','r'); hold on;
plot(t, LBMPC50.sysH(4,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time ');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(t, LMPC50.sysH(5,:),'Linewidth',1.5,'Color','r'); hold on;
plot(t, LBMPC50.sysH(5,:),'Linewidth',1.5,'Color','b'); hold on;
plot(t,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('time ');
ylabel('\delta u');
title('Sys input');


figure;
plot(LMPC50.sysH(1,:),LMPC50.sysH(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
plot(LBMPC50.sysH(1,:),LBMPC50.sysH(2,:),'Linewidth',1.5,'Marker','.','Color','b'); 
grid on
xlabel('\delta x_1');
ylabel('\delta x_2');
title('State space');

legend({'LMPC N=50','LBMPC N=50'},'Location','northwest')
