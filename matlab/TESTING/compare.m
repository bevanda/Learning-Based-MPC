%% Load data
LMPC_sys=load('LMPC_sysH_K_T01_N40_new_new.mat');
LMPC_art=load('LMPC_artH_K_T01_N40_new_new.mat');
% NMPC_sys=load('NMPC_sysH_K_T01_N20_new_1000iter.mat');
% NMPC_art=load('NMPC_artH_K_T01_N20_new_1000iter.mat');
NMPC_sys=load('NMPC_sysH_K_T01_N40_new.mat');
NMPC_art=load('NMPC_artH_K_T01_N40_new.mat');
% NMPC_sys=load('LMPC_sysH_K_T01_N20_new_1000iter_T1000P.mat');
% NMPC_art=load('LMPC_artH_K_T01_N20_new_1000iter_T1000P.mat');
%% Plot
iterations=size(LMPC_sys.sysHistory,2)-1;
iterations2=size(NMPC_sys.sysHistory,2)-1;
figure;
subplot(5,1,1);
plot(0:iterations2,NMPC_sys.sysHistory(1,:), 0:iterations,LMPC_sys.sysHistory(1,:), 'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('x1');

title('mass flow');
legend({ 'NMPC', 'LMPC'},'Location','northwest')
subplot(5,1,2);
plot(0:iterations2,NMPC_sys.sysHistory(2,:),0:iterations,LMPC_sys.sysHistory(2,:), 'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('x2');
title('pressure rise');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')
subplot(5,1,3);
plot( 0:iterations2,NMPC_sys.sysHistory(3,:),0:iterations,LMPC_sys.sysHistory(3,:),'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('x3');
title('throttle');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')
subplot(5,1,4);
plot( 0:iterations2,NMPC_sys.sysHistory(4,:),0:iterations,LMPC_sys.sysHistory(4,:),'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('x4');
title('throttle rate');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')
subplot(5,1,5);
plot( 0:iterations2,NMPC_sys.sysHistory(5,:), 0:iterations,LMPC_sys.sysHistory(5,:),'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('u');
title('input');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')

figure;
plot_refs=plot(0:iterations,LMPC_art.art_refHistory(1,:), 0:iterations,LMPC_sys.sysHistory(1,:),'Linewidth',2); hold on;
plot_refs(1).LineStyle='--';
plot_refs(1).Color='m';
plot_refs(2).Color='r';
plot_refs2=plot(0:iterations2,NMPC_art.art_refHistory(1,:), 0:iterations2,NMPC_sys.sysHistory(1,:),'Linewidth',2);
plot_refs2(1).LineStyle='-.';
plot_refs2(1).Color='c';
plot_refs2(2).Color='b';
grid on
xlabel('iterations');
% ylabel('references');
title('Mass flow relevant responses');
legend({'art ref LMPC', 'mass flow LMPC','art ref NMPC', 'mass flow NMPC'},'Location','southeast')


figure;
plot(NMPC_sys.sysHistory(1,:),NMPC_sys.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','b'); hold on;
plot(LMPC_sys.sysHistory(1,:),LMPC_sys.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r');
grid on
xlabel('x1');
ylabel('x2');
title('State space');
legend({'NMPC', 'LMPC'},'Location','southeast')
