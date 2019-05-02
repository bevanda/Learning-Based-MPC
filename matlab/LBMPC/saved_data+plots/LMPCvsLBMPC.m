close all;
addpath('./data/')
%% Load data
LMPC50=load('LMPC_N50_sys_k1cons.mat');
LBMPC50=load('LBMPC_N50_sys_newP_newpoly_k1cons.mat');
LBMPC60=load('LBMPC_N60_sys_newP_newpoly_k1cons.mat');

%% Plot
iterations=size(LBMPC50.sysH,2)-1;

figure;
subplot(5,1,1);
plot(0:iterations, LMPC50.sysH(1,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC50.sysH(1,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations,LBMPC60.sysH(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');

subplot(5,1,2);
plot(0:iterations, LMPC50.sysH(2,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC50.sysH(2,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations,LBMPC60.sysH(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');
legend({'LMPC N=50','LBMPC N=50', 'LBMPC N=60'},'Location','northwest')

subplot(5,1,3);
plot(0:iterations, LMPC50.sysH(3,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC50.sysH(3,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations,LBMPC60.sysH(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(0:iterations, LMPC50.sysH(4,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC50.sysH(4,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations,LBMPC60.sysH(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(0:iterations, LMPC50.sysH(5,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC50.sysH(5,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations,LBMPC60.sysH(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');


figure;
plot(LMPC50.sysH(1,:),LMPC50.sysH(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
plot(LBMPC50.sysH(1,:),LBMPC50.sysH(2,:),'Linewidth',1.5,'Marker','.','Color','b'); hold on;
plot(LBMPC60.sysH(1,:),LBMPC60.sysH(2,:),'Linewidth',1.5,'Marker','.','Color','g'); 

grid on
xlabel('\delta x_1');
ylabel('\delta x_2');
title('State space');

legend({'LMPC N=50','LBMPC N=50', 'LBMPC N=60'},'Location','northwest')
