%% Load data
LMPC_sys=load('LMPC_N50_sys.mat');
% LMPC_art=load('LMPC_artH_K_T01_N40_new_new.mat');
% NMPC_sys=load('NMPC_sysH_K_T01_N20_new_1000iter.mat');
% NMPC_art=load('NMPC_artH_K_T01_N20_new_1000iter.mat');
LBMPC_sys=load('LBMPC_N60_sys.mat');
% NMPC_art=load('NMPC_artH_K_T01_N40_new.mat');
% LBMPC_sys=load('LMPC_sysH_K_T01_N20_new_1000iter_T1000P.mat');
% NMPC_art=load('LMPC_artH_K_T01_N20_new_1000iter_T1000P.mat');
%% Plot
iterations=size(LMPC_sys.sysHistory,2)-1;
iterations2=size(LBMPC_sys.sysHistory,2)-1;
figure;
subplot(5,1,1);
plot(0:iterations2,LBMPC_sys.sysHistory(1,:), 0:iterations,LMPC_sys.sysHistory(1,:), 'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('\delta x_1');

title('mass flow');
legend({ 'LBMPC', 'LMPC'},'Location','northwest')
subplot(5,1,2);
plot(0:iterations2,LBMPC_sys.sysHistory(2,:),0:iterations,LMPC_sys.sysHistory(2,:), 'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('\delta x_2');
title('pressure rise');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')
subplot(5,1,3);
plot( 0:iterations2,LBMPC_sys.sysHistory(3,:),0:iterations,LMPC_sys.sysHistory(3,:),'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('\delta x_3');
title('throttle');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')
subplot(5,1,4);
plot( 0:iterations2,LBMPC_sys.sysHistory(4,:),0:iterations,LMPC_sys.sysHistory(4,:),'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('\delta x_4');
title('throttle rate');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')
subplot(5,1,5);
plot( 0:iterations2,LBMPC_sys.sysHistory(5,:), 0:iterations,LMPC_sys.sysHistory(5,:),'Linewidth',1.5);
grid on
xlabel('iterations');
ylabel('\delta u');
title('input');
% legend({ 'NMPC', 'LMPC'},'Location','southeast')

% figure;
% plot_refs=plot(0:iterations,LMPC_art.art_refHistory(1,:), 0:iterations,LMPC_sys.sysHistory(1,:),'Linewidth',2); hold on;
% plot_refs(1).LineStyle='--';
% plot_refs(1).Color='m';
% plot_refs(2).Color='r';
% plot_refs2=plot(0:iterations2,NMPC_art.art_refHistory(1,:), 0:iterations2,LBMPC_sys.sysHistory(1,:),'Linewidth',2);
% plot_refs2(1).LineStyle='-.';
% plot_refs2(1).Color='c';
% plot_refs2(2).Color='b';
% grid on
% xlabel('iterations');
% % ylabel('references');
% title('Mass flow relevant responses');
% legend({'art ref LMPC', 'mass flow LMPC','art ref NMPC', 'mass flow NMPC'},'Location','southeast')


figure;
plot(LBMPC_sys.sysHistory(1,:),LBMPC_sys.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','b'); hold on;
plot(LMPC_sys.sysHistory(1,:),LMPC_sys.sysHistory(2,:),'Linewidth',1.5,'Marker','.','Color','r');
grid on
xlabel('x1');
ylabel('x2');
title('State space');
legend({'LBMPC', 'LMPC'},'Location','southeast')
%% Qualitative euvaluation of the responses
disp('OVERSHOOT\n');
disp('Evaluating LBMPC responses\n');
fprintf('Percent of mass flow over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(1,:),0));
fprintf('Percent of pressure rise over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(2,:),0));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(4,:),0));
disp('Evaluating LMPC responses\n');
fprintf('Percent of mass flow over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(1,:),0));
fprintf('Percent of pressure rise over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(2,:),0));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(4,:),0));
%%
disp('SETTLING TIME\n');
disp('Evaluating LBMPC responses\n');
t_ss_lbmpc1=compute_ss_time(0.01*(0:iterations),LBMPC_sys.sysHistory(1,:),0,0.02);
fprintf('Settling time mass flow: %d\n',t_ss_lbmpc1);
t_ss_lbmpc2=compute_ss_time(0.01*(0:iterations),LBMPC_sys.sysHistory(2,:),0,0.02);
fprintf('Settling time pressure rise: %d\n',t_ss_lbmpc2);
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(4,:),0));
disp('Evaluating LMPC responses\n');
t_ss_lmpc1=compute_ss_time(0.01*(0:iterations),LMPC_sys.sysHistory(1,:),0,0.02);
fprintf('Settling time mass flow: %d\n',t_ss_lmpc1);
t_ss_lmpc2=compute_ss_time(0.01*(0:iterations),LMPC_sys.sysHistory(2,:),0,0.02);
fprintf('Settling time pressure rise: %d\n',t_ss_lmpc2);
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(4,:),0));
%%
disp('CTRL ERROR\n');
disp('Evaluating LBMPC responses\n');
fprintf('Accmumulated ctrl error mass flow: %d\n',compute_ctrl_error(0.01*(0:iterations),LBMPC_sys.sysHistory(1,:),0,t_ss_lbmpc1));
fprintf('Accmumulated ctrl error pressure rise: %d\n',compute_ctrl_error(0.01*(0:iterations),LBMPC_sys.sysHistory(2,:),0,t_ss_lbmpc2));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(4,:),0));
disp('Evaluating LMPC responses\n');
fprintf('Accmumulated ctrl error mass flow: %d\n',compute_ctrl_error(0.01*(0:iterations),LMPC_sys.sysHistory(1,:),0,t_ss_lmpc1));
fprintf('Accmumulated ctrl error pressure rise: %d\n',compute_ctrl_error(0.01*(0:iterations),LMPC_sys.sysHistory(2,:),0,t_ss_lmpc2));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(4,:),0));
%%
disp('CTRL energy\n');
disp('Evaluating LBMPC responses\n');
fprintf('Accmumulated ctrl effort mass flow: %d\n',compute_ctrl_energy(0.01*(0:iterations),LBMPC_sys.sysHistory(5,:),t_ss_lbmpc1));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(4,:),0));
disp('Evaluating LMPC responses\n');
fprintf('Accmumulated ctrl effort mass flow: %d\n',compute_ctrl_energy(0.01*(0:iterations),LMPC_sys.sysHistory(5,:),t_ss_lmpc1));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(4,:),0));

%%
eps=0.05;
disp('setl time\n');
disp('Evaluating LBMPC responses\n');
fprintf('setl time mass flow: %d\n',compute_ss_time(0.01*(0:iterations),LBMPC_sys.sysHistory(1,:),0,eps));
fprintf('setl time pressure rise: %d\n',compute_ss_time(0.01*(0:iterations),LBMPC_sys.sysHistory(2,:),0,eps));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(4,:),0));
disp('Evaluating LMPC responses\n');
fprintf('Accmumulated ctrl error mass flow: %d\n',compute_ss_time(0.01*(0:iterations),LMPC_sys.sysHistory(1,:),0,eps));
fprintf('Accmumulated ctrl error pressure rise: %d\n',compute_ss_time(0.01*(0:iterations),LMPC_sys.sysHistory(2,:),0,eps));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(4,:),0));