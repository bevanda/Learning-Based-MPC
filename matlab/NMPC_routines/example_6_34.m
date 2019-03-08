function example_6_34
% example_6_34
% For this example, Assumption 6.4 does not hold. As a consequence,
% Theorem 6.21 is not applicalbe and we cannot expect asymptotic stability
% of the closed loop. As Theorem 6.33 predicts, the solutions converge to
% smaller and smaller neighbourhoods of the target as N increases, but they
% do not converge to the target for fixed N.
    addpath('./nmpcroutine');
    clear all;
    close all;

    mpciterations = 10;
    N             = 2;
    T             = 0.1;
    tmeasure      = 0.0;
    xmeasure      = [2.0];
    u0            = zeros(1,N);
    iprint        = 5;
    tol_opt       = 1e-8;
    opt_option    = 0;
    type          = 'difference equation';
    atol_ode_real = [];
    rtol_ode_real = [];
    atol_ode_sim  = [];
    rtol_ode_sim  = [];

    nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system, ...
         mpciterations, N, T, tmeasure, xmeasure, u0, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint, @printHeader, @printClosedloopData, @plotTrajectories);

    rmpath('./nmpcroutine');
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)
    cost = x^2+abs(u);
end

function cost = terminalcosts(t, x)
    cost = 0.0;
end

function [c,ceq] = constraints(t, x, u)
    c   = [];
    ceq = [];
end

function [c,ceq] = terminalconstraints(t, x)
    c   = [];
    ceq = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = [];
    ub  = [];
end

function y = system(t, x, u, T)
    y(1) = x(1)+u(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of output format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printHeader()
    fprintf('   k  |      u(k)        x(1)        x(2)     Time\n');
    fprintf('--------------------------------------------------\n');
end

function printClosedloopData(mpciter, u, x, t_Elapsed)
    fprintf(' %3d  | %+11.6f %+11.6f  %+6.3f', mpciter, u, x, t_Elapsed);
end

function plotTrajectories(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
    figure(1);
        title('Closed loop trajectory');
        xlabel('t');
        ylabel('x');
        grid on;
        hold on;
        plot(t_intermediate(:),x_intermediate(:),'-ok');
        axis([0 1 0 2]);
        axis square;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
