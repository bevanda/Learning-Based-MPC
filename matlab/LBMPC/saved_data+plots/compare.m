%% Load data
LBMPC_N80_sys=load('LBMPC_N80_sys.mat');
LBMPC_N100_sys=load('LBMPC_N100_sys.mat');
%% Qualitative assesment of responses


%%
disp('CTRL ERROR\n');
disp('Evaluating LBMPC 80responses\n');
fprintf('Accmumulated ctrl error mass flow: %d\n',compute_ctrl_error(0.01*(0:iterations),LBMPC_N100_sys.sysHistory(1,:),0,t_ss_80lbmpc1));
fprintf('Accmumulated ctrl error pressure rise: %d\n',compute_ctrl_error(0.01*(0:iterations),LBMPC_N100_sys.sysHistory(2,:),0,t_ss_80lbmpc2));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(4,:),0));
disp('Evaluating LBMPC 100 responses\n');
fprintf('Accmumulated ctrl error mass flow: %d\n',compute_ctrl_error(0.01*(0:iterations),LBMPC_N80_sys.sysHistory(1,:),0,t_ss_100lmpc1));
fprintf('Accmumulated ctrl error pressure rise: %d\n',compute_ctrl_error(0.01*(0:iterations),LBMPC_N80_sys.sysHistory(2,:),0,t_ss_100lmpc2));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(4,:),0));
%%
disp('CTRL energy\n');
disp('Evaluating LBMPC 80 responses\n');
fprintf('Accmumulated ctrl effort mass flow: %d\n',compute_ctrl_energy(0.01*(0:iterations),LBMPC_N100_sys.sysHistory(5,:),t_ss_lbmpc1));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LBMPC_sys.sysHistory(4,:),0));
disp('Evaluating LBMPC 100 responses\n');
fprintf('Accmumulated ctrl effort mass flow: %d\n',compute_ctrl_energy(0.01*(0:iterations),LBMPC_N80_sys.sysHistory(5,:),t_ss_lmpc1));
% fprintf('Percent of throttle over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(3,:),0));
% fprintf('Percent of throttle rate over/under-shoot: %d\n',compute_overshoot(LMPC_sys.sysHistory(4,:),0));
