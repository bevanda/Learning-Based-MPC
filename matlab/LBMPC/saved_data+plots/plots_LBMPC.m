close all;
addpath('./data/')
%% Load data
LMPC50=load('LMPC_N50_sys.mat');
LBMPC60=load('LBMPC_N60_sys_newP_newpoly.mat');
LBMPC80=load('LBMPC_N80_sys_newP.mat');
LBMPC100=load('LBMPC_N100_sys_newP.mat');

%% Plot
iterations=size(LBMPC60.sysH,2)-1;

figure;
subplot(5,1,1);
plot(0:iterations, LMPC50.sysHistory(1,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC60.sysH(1,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80.sysHistory(1,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LBMPC100.sysHistory(1,:),'Linewidth',1.5,'Color','b'); hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_1');
title('mass flow');
legend({'LMPC','LBMPC60', 'LBMPC N=80', 'LBMPC N=100'},'Location','northwest')

subplot(5,1,2);
plot(0:iterations, LMPC50.sysHistory(2,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC60.sysH(2,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80.sysHistory(2,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LBMPC100.sysHistory(2,:),'Linewidth',1.5,'Color','b');  hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');

subplot(5,1,3);
plot(0:iterations, LMPC50.sysHistory(3,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC60.sysH(3,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80.sysHistory(3,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LBMPC100.sysHistory(3,:),'Linewidth',1.5,'Color','b');  hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');

subplot(5,1,4);
plot(0:iterations, LMPC50.sysHistory(4,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC60.sysH(4,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80.sysHistory(4,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LBMPC100.sysHistory(4,:),'Linewidth',1.5,'Color','b');  hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');

subplot(5,1,5);
plot(0:iterations, LMPC50.sysHistory(5,:),'Linewidth',1.5,'Color','r'); hold on;
plot(0:iterations, LBMPC60.sysH(5,:),'Linewidth',1.5,'Color','m'); hold on;
plot(0:iterations,LBMPC80.sysHistory(5,:),'Linewidth',1.5,'Color','g'); hold on;
plot(0:iterations,LBMPC100.sysHistory(5,:),'Linewidth',1.5,'Color','b');  hold on;
plot(0:iterations,0*(0:iterations),'Linewidth',1,'LineStyle',':','Color','k'); 
grid on;
xlabel('iterations');
ylabel('\delta u');
title('Sys input');


figure;
plot(LMPC50.sysHistory(1,:),LMPC50.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r'); hold on;
plot(LBMPC60.sysH(1,:),LBMPC60.sysH(2,:),'Linewidth',1.5,'Marker','.','Color','m'); hold on;
plot(LBMPC80.sysHistory(1,:),LBMPC80.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','g'); hold on;
plot(LBMPC100.sysHistory(1,:),LBMPC100.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','b'); 

grid on
xlabel('\delta x_1');
ylabel('\delta x_2');
title('State space');

legend({'LMPC', 'LBMPC60','LBMPC N=80', 'LBMPC N=100'},'Location','northwest')
