function example_fig_7_1
% example_fig_7_1
% This example shows an evolution of "alpha" values if the weight of the
% terminal cost function term is increased successively.
    addpath('./mpcalpha');
    clearvars;
    close all;

    C         = 2;
    sigma     = 0.55;
    N         = 5;
    m         = 1;
    endweight = 1:20;

    alpha_omega_plot(C, sigma, N, m, endweight);
    
    rmpath('./mpcalpha');
end

