close all;
addpath('./data/')
%% Load data
LBMPC50r=load('LBMPC_N50_sys_full');
LBMPC60r=load('LBMPC_N60_sys_full');
LBMPC80r=load('LBMPC_N80_sys_full');

%% Plot
iterations=size(LBMPC80r.sysH,2)-1;

figure;
subplot(5,1,1);
plot(0:iterations, LBMPC50r.sysH(1,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LBMPC60r.sysH(1,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80r.sysH(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');


subplot(5,1,2);
plot(0:iterations, LBMPC50r.sysH(2,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LBMPC60r.sysH(2,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations, LBMPC80r.sysH(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');
legend({'LBMPC N=50' 'LBMPC N=60','LBMPC N=80'},'Location','northwest')

subplot(5,1,3);
plot(0:iterations, LBMPC50r.sysH(3,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LBMPC60r.sysH(3,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80r.sysH(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(0:iterations, LBMPC50r.sysH(4,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LBMPC60r.sysH(4,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80r.sysH(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(0:iterations, LBMPC50r.sysH(5,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations, LBMPC60r.sysH(5,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80r.sysH(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');



