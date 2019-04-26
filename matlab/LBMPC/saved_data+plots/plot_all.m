close all;
addpath('./data/');
%% Load data
LBMPC40=load('LBMPC_N40_sys_CORRECT.mat');
NMPC40=load('NMPC_N40_sys_CORRECT.mat');
LMPC40=load('LMPC_N40_sys.mat');
%% Plot
iterations=size(LBMPC40.sysHistory,2)-1;

figure;
subplot(5,1,1);
plot(0:iterations, LBMPC40.sysHistory(1,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,NMPC40.sysHistory(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC40.sysHistory(1,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');
legend({'LBMPC','NMPC', 'LMPC'},'Location','northwest')

subplot(5,1,2);
plot(0:iterations, LBMPC40.sysHistory(2,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,NMPC40.sysHistory(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC40.sysHistory(2,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');

subplot(5,1,3);
plot(0:iterations, LBMPC40.sysHistory(3,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,NMPC40.sysHistory(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC40.sysHistory(3,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');


subplot(5,1,4);
plot(0:iterations, LBMPC40.sysHistory(4,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,NMPC40.sysHistory(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC40.sysHistory(4,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');


subplot(5,1,5);
plot(0:iterations, LBMPC40.sysHistory(5,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,NMPC40.sysHistory(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC40.sysHistory(5,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');



figure;
plot(LBMPC40.sysHistory(1,:),LBMPC40.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
plot(NMPC40.sysHistory(1,:),NMPC40.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','g'); hold on;
plot(LMPC40.sysHistory(1,:),LMPC40.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','b'); 

grid on
xlabel('\delta x_1');
ylabel('\delta x_2');
title('State space');

legend({'LBMPC','NMPC', 'LMPC'},'Location','northwest')
