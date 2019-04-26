%% Load data
close all;
addpath('./data/');

LBMPC100_q100=load('LBMPC_N100_sys.mat');
LBMPC100_q200=load('LBMPC_N100_sys_q200.mat');
LBMPC100_q500=load('LBMPC_N100_sys_q500.mat');
NMPC100=load('NMPC_N100_sys.mat');
%% Plot
iterations=size(NMPC100.sysHistory,2)-1;
figure;
subplot(5,1,1);
plot(0:iterations,LBMPC100_q100.sysHistory(1,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q200.sysHistory(1,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q500.sysHistory(1,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,NMPC100.sysHistory(1,:),'Linewidth',1.5,'LineStyle','-');
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');
legend({'LBMPC q=100', 'LBMPC q=200', 'LBMPC q=500', 'NMPC'},'Location','northwest')

subplot(5,1,2);
plot(0:iterations,LBMPC100_q100.sysHistory(2,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q200.sysHistory(2,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q500.sysHistory(2,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,NMPC100.sysHistory(2,:),'Linewidth',1.5,'LineStyle','-');
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');

subplot(5,1,3);
plot(0:iterations,LBMPC100_q100.sysHistory(3,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q200.sysHistory(3,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q500.sysHistory(3,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,NMPC100.sysHistory(3,:),'Linewidth',1.5,'LineStyle','-');
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(0:iterations,LBMPC100_q100.sysHistory(4,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q200.sysHistory(4,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q500.sysHistory(4,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,NMPC100.sysHistory(4,:),'Linewidth',1.5,'LineStyle','-');
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(0:iterations,LBMPC100_q100.sysHistory(5,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q200.sysHistory(5,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,LBMPC100_q500.sysHistory(5,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
plot(0:iterations,NMPC100.sysHistory(5,:),'Linewidth',1.5,'LineStyle','-');
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');

