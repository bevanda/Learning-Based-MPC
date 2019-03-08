function example_6_20
% example_6_20
% This example shows the improvement in "alpha" for unchanged values of "C"
% and "sigma" if the horizon is increased.
    addpath('./mpcalpha');

    C = 1;
    sigma = exp(-1);
    for N = 2:3;
        alpha = mpcalpha_Csigma(C,sigma,N);
        fprintf(['Horizon: %3d   alpha = %+10.6f   Tradeoff: %+10.6f', ...
                 '\n'], N, alpha, (100/alpha)-100);
    end

    rmpath('./mpcalpha');
end

