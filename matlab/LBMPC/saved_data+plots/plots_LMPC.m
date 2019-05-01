close all;
addpath('./data/')
%% Load data
LMPC=load('LMPC_N20_sys_full.mat');
LMPCrk1=load('LMPC_N20_sys_k1cons.mat');
LMPCrT=load('LMPC_N20_sys.mat');


%% Plot
iterations=size(LMPC.sysH,2)-1;

figure;
subplot(5,1,1);
plot(0:iterations, LMPC.sysH(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations, LMPCrk1.sysH(1,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LMPCrT.sysH(1,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');


subplot(5,1,2);
plot(0:iterations, LMPC.sysH(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations, LMPCrk1.sysH(2,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LMPCrT.sysH(2,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');
legend({'LMPC' 'LMPCr k1','LMPCr N'},'Location','northwest')

subplot(5,1,3);
plot(0:iterations, LMPC.sysH(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations, LMPCrk1.sysH(3,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LMPCrT.sysH(3,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(0:iterations, LMPC.sysH(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations, LMPCrk1.sysH(4,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LMPCrT.sysH(4,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(0:iterations, LMPC.sysH(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations, LMPCrk1.sysH(5,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LMPCrT.sysH(5,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');



