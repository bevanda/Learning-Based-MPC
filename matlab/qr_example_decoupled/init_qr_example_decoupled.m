%==========================================================================
% Configuration for MPT toolbox
%==========================================================================
clc
close all
clear all

if strcmp(getenv('USER'), 'bouffard') == 0
    % Stuff for Anil only
    mpt_path = 'C:\Program Files\MATLAB\R2009b\toolbox\mpt';
    disp('Adding MPT toolbox to MATLAB path...');
    addpath(genpath(mpt_path));
    
    %% Misc configuration
    quad_dat_fname = 'quad.dat';
    quad_dat_fname = 'quad.mat';
    dual_ekf_fname = 'dual_ekf.dat';
    
    % added by George
    
	% quad_bin_fname = 'quad.bin';
else
    % Stuff for Pat only
    clear all;
    %clc;
    close all;
    
    %% Misc configuration
    quad_bin_fname = 'quad_xy.bin';
end

N = 15; % MPC horizon

quad_bin_fname = 'quad_xy.bin';
define_system_xy;
init2; % most calcs are done here
init3; % write out to file(s)
