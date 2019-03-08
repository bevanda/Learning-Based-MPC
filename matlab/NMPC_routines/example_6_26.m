function example_6_26
% example_6_26
% This example show the impact of the horizon length on the stability of
% the closed loop. Here, the closed loop converges towards the desired
% equilibrium (0, -0.5) for N > 10, but the solution remains unchanged
% for N <= 10.
    addpath('./nmpcroutine');
    clear all;
    close all;

    mpciterations = 25;
    N             = 11;
    T             = 0.1;
    tmeasure      = 0.0;
    xmeasure      = [0.0, 0.5];
    u0            = 0.2*ones(1,N);
    tol_opt       = 1e-8;
    opt_option    = 0;
    iprint        = 5;
    type          = 'difference equation';
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;

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
    cost = norm(x-[0, -0.5],2)^2+norm(u,2)^2;
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
    lb  = 0;
    ub  = 0.2;
end

function y = system(t, x, u, T)
    xx = [x(1), 2*x(2)];
    n = norm(xx,2);

    if (n == 0.0)
        dx = [0.0 0.0];
    else
        phi = real(acos(xx(2)/n));
        if (xx(1) < 0)
        phi = 2*pi-phi;
        end

        y(1) = n*sin(phi+u(1));
        y(2) = n*cos(phi+u(1))/2;
    end
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
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f  %+6.3f', mpciter, ...
            u(1,1), x(1), x(2), t_Elapsed);
end

function plotTrajectories(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
    figure(1);
        title('x_1/x_2 closed loop trajectory');
        xlabel('x_1');
        ylabel('x_2');
        grid on;
        hold on;
        plot(sin(0:pi/20:pi), cos(0:pi/20:pi)/2, ...
             'Color',[0.8 0.8 0.8],'LineWidth',2);
        plot(0, -0.5,'ok','MarkerSize',8);
        plot(x_intermediate(:,1),x_intermediate(:,2),'or', ...
             'MarkerFaceColor','r');
        axis([-0.5 1.5 -1 1]);
        axis square;

    figure(2);
        title(['x_1 and x_2 closed loop trajectory']);
        xlabel('n');
        ylabel('x_1(n), x_2(n)');
        grid on;
        hold on;
        plot(t_intermediate,x_intermediate(:,1),'-ok');
        plot(t_intermediate,x_intermediate(:,2),'-ok');
        axis([0 2.5 -0.5 1]);
        axis square;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
