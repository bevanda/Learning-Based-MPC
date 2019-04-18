%% Load data
close all;
LBMPC60=load('LBMPC_N60_sys.mat');
LBMPC40=load('LBMPC_N50_sys.mat');
LMPC=load('LMPC_N50_sys.mat');

%% Plot
iterations=size(LMPC.sysHistory,2)-1;

figure;
subplot(5,1,1);
plot(0:iterations, LBMPC60.sysHistory(1,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC40.sysHistory(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(1,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');
legend({ 'LBMPC N=60', 'LBMPC N=50', 'LMPC N=50'},'Location','northwest')

subplot(5,1,2);
plot(0:iterations, LBMPC60.sysHistory(2,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC40.sysHistory(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(2,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');
% legend({ 'LBMPC40','LBMPC40', 'LMPC'},'Location','southeast')

subplot(5,1,3);
plot(0:iterations, LBMPC60.sysHistory(3,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC40.sysHistory(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(3,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');
% legend({ 'LBMPC40', 'LBMPC40','LMPC'},'Location','southeast')

subplot(5,1,4);
plot(0:iterations, LBMPC60.sysHistory(4,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC40.sysHistory(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(4,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');
% legend({ 'LBMPC40','LBMPC40', 'LMPC'},'Location','southeast')

subplot(5,1,5);
plot(0:iterations, LBMPC60.sysHistory(5,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,LBMPC40.sysHistory(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LMPC.sysHistory(5,:),'Linewidth',1.5,'Color','b'); 
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');
% legend({ 'LBMPC40','LBMPC40', 'LMPC'},'Location','southeast')


figure;
plot(LBMPC60.sysHistory(1,:),LBMPC60.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
plot(LBMPC40.sysHistory(1,:),LBMPC40.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','g'); hold on;
plot(LMPC.sysHistory(1,:),LMPC.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','b'); 

grid on
xlabel('\delta x_1');
ylabel('\delta x_2');
title('State space');

legend({'LBMPC60','LBMPC40', 'LMPC'},'Location','southeast')
