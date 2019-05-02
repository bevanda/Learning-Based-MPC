%% Load data
close all;
% LBMPC60=load('LBMPC_N60_sys.mat');
LBMPC50=load('LBMPC_N50_sys_newP_newpoly.mat');
LBMPC40=load('LBMPC_N40_sys_CORRECT.mat');
% NMPC60=load('NMPC_N60_sys.mat');
NMPC50=load('NMPC_N50_sys_newP.mat');
NMPC40=load('NMPC_N40_sys_CORRECT.mat');
%% Plot
iterations=size(LBMPC50.sysH,2)-1;

figure;
subplot(5,1,1);
% plot(0:iterations, LBMPC60.sysHistory(1,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC50.sysH(1,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC40.sysHistory(1,:),'Linewidth',1.5,'Color','b','LineStyle','-'); hold on;
% plot(0:iterations, NMPC60.sysHistory(1,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC50.sysH(1,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC40.sysHistory(1,:),'Linewidth',1.5,'Color','b','LineStyle','-.');
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');
legend({'LBMPC N=80', 'LBMPC N=100', 'NMPC N=80', 'NMPC N=100'},'Location','northwest')

subplot(5,1,2);
% plot(0:iterations, LBMPC60.sysHistory(2,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC50.sysH(2,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC40.sysHistory(2,:),'Linewidth',1.5,'Color','b','LineStyle','-'); hold on;
% plot(0:iterations, NMPC60.sysHistory(2,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC50.sysH(2,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC40.sysHistory(2,:),'Linewidth',1.5,'Color','b','LineStyle','-.'); 
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');
% legend({ 'LBMPC40','LBMPC40', 'LMPC'},'Location','southeast')

subplot(5,1,3);
% plot(0:iterations, LBMPC60.sysHistory(3,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC50.sysH(3,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC40.sysHistory(3,:),'Linewidth',1.5,'Color','b','LineStyle','-'); hold on;
% plot(0:iterations, NMPC60.sysHistory(3,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC50.sysH(3,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC40.sysHistory(3,:),'Linewidth',1.5,'Color','b','LineStyle','-.');
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');
% legend({ 'LBMPC40', 'LBMPC40','LMPC'},'Location','southeast')

subplot(5,1,4);
% plot(0:iterations, LBMPC60.sysHistory(4,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC50.sysH(4,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC40.sysHistory(4,:),'Linewidth',1.5,'Color','b','LineStyle','-'); hold on;
% plot(0:iterations, NMPC60.sysHistory(4,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC50.sysH(4,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC40.sysHistory(4,:),'Linewidth',1.5,'Color','b','LineStyle','-.');
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');
% legend({ 'LBMPC40','LBMPC40', 'LMPC'},'Location','southeast')

subplot(5,1,5);
% plot(0:iterations, LBMPC60.sysHistory(5,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC50.sysH(5,:),'Linewidth',1.5,'Color','r','LineStyle','-'); hold on;
plot(0:iterations,LBMPC40.sysHistory(5,:),'Linewidth',1.5,'Color','b','LineStyle','-'); hold on;
% plot(0:iterations, NMPC60.sysHistory(5,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC50.sysH(5,:),'Linewidth',1.5,'Color','r','LineStyle','-.'); hold on;
plot(0:iterations,NMPC40.sysHistory(5,:),'Linewidth',1.5,'Color','b','LineStyle','-.');
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');
% legend({ 'LBMPC40','LBMPC40', 'LMPC'},'Location','southeast')


% figure;
% plot(LBMPC60.sysHistory(1,:),LBMPC60.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
% plot(LBMPC80.sysHistory(1,:),LBMPC80.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
% plot(LBMPC100.sysHistory(1,:),LBMPC100.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','b'); 
% 
% grid on
% xlabel('\delta x_1');
% ylabel('\delta x_2');
% title('State space');

% legend({'LBMPC60','LBMPC80', 'LBMPC100'},'Location','southeast')
