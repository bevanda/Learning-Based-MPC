%% Load data
close all;
NMPC=load('NMPC_sysH_K_T01_N40_new.mat');
LBMPC=load('LBMPC_N40_sys.mat');
LMPC=load('LMPC_N40_sys.mat');

%% Plot
iterations=size(LMPC.sysHistory,2)-1;

figure;
subplot(5,1,1);
plot(0:iterations, NMPC.sysHistory(1,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC.sysHistory(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(1,:),'Linewidth',1.5,'Color','b'); 
grid on
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');
legend({ 'NMPC', 'LBMPC', 'LMPC'},'Location','northwest')

subplot(5,1,2);
plot(0:iterations, NMPC.sysHistory(2,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC.sysHistory(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(2,:),'Linewidth',1.5,'Color','b'); 
grid on
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');
% legend({ 'NMPC','LBMPC', 'LMPC'},'Location','southeast')

subplot(5,1,3);
plot(0:iterations, NMPC.sysHistory(3,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC.sysHistory(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(3,:),'Linewidth',1.5,'Color','b'); 
grid on
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');
% legend({ 'NMPC', 'LBMPC','LMPC'},'Location','southeast')

subplot(5,1,4);
plot(0:iterations, NMPC.sysHistory(4,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC.sysHistory(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(4,:),'Linewidth',1.5,'Color','b'); 
grid on
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');
% legend({ 'NMPC','LBMPC', 'LMPC'},'Location','southeast')

subplot(5,1,5);
plot(0:iterations, NMPC.sysHistory(5,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC.sysHistory(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(5,:),'Linewidth',1.5,'Color','b'); 
grid on
xlabel('iterations');
ylabel('\delta u');
title('Sys input');
% legend({ 'NMPC','LBMPC', 'LMPC'},'Location','southeast')




figure;
plot(NMPC.sysHistory(1,:),NMPC.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
plot(LBMPC.sysHistory(1,:),LBMPC.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','g'); hold on;
plot(LMPC.sysHistory(1,:),LMPC.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','b'); 

grid on
xlabel('\delta x_1');
ylabel('\delta x_2');
title('State space');
legend({'NMPC','LBMPC', 'LMPC'},'Location','southeast')


