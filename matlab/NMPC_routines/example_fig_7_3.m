function example_fig_7_3
% example_fig_7_3
% This example shows an evolution of "alpha" values if the control horizon
% is increased, i.e. if successively all the open loop controls are
% implemented.
    addpath('./mpcalpha');
    clearvars;
    close all;

    C         = 2;
    sigma     = 0.75;
    N         = 11;
    m         = 1:10;
    
    alpha_m_plot(C, sigma, N, m);
    
    rmpath('./mpcalpha');
end

