%% Qualitative assesment of responses
% load data
clc; close all;
addpath('./data/');
addpath('./utilities/');

LBMPC_N50_sys=load('LBMPC_N50_sys_full.mat');
LBMPC_N60_sys=load('LBMPC_N60_sys_full.mat');
LBMPC_N80_sys=load('LBMPC_N80_sys_full.mat');
Ts=0.01;
iterations=length(LBMPC_N50_sys.sysH)-1;
t=Ts*(0:iterations);

% get response info
disp('LBMPC N=50 performance');
disp('x_1');
response_info(LBMPC_N50_sys.sysH(1,:),LBMPC_N50_sys.sysH(5,:),t,0)
disp('x_2');
response_info(LBMPC_N50_sys.sysH(2,:),LBMPC_N50_sys.sysH(5,:),t,0)

disp('LBMPC N=60 performance');
disp('x_1');
response_info(LBMPC_N60_sys.sysH(1,:),LBMPC_N60_sys.sysH(5,:),t,0)
disp('x_2');
response_info(LBMPC_N60_sys.sysH(2,:),LBMPC_N60_sys.sysH(5,:),t,0)

disp('LBMPC N=80 performance');
disp('x_1');
response_info(LBMPC_N80_sys.sysH(1,:),LBMPC_N80_sys.sysH(5,:),t,0)
disp('x_2');
response_info(LBMPC_N80_sys.sysH(2,:),LBMPC_N80_sys.sysH(5,:),t,0)