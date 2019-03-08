function example_8_31
% example_8_31
% This example implements Artstein's circles. The boundary for the state
% solution can be shifted by setting the variable bound. Moreover, both the
% continuous time and the discrete time version are implemented and can be
% switched by supplying "system_ct" instead of "system_dt" to the NMPC routine.
% Note that the variable "type" needs to be reset as well.
    addpath('./nmpcroutine');
    clear all;
    close all;

    global bound;

    bound         = 0.1;

    mpciterations = 30;
    N             = 2;
    T             = 1;
    tmeasure      = 0.0;
    xmeasure      = [0, 1];
    u0            = ones(1,N);
    tol_opt       = 1e-8;
    opt_option    = 0;
    iprint        = 5;
    type          = 'difference equation';
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;

    nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system_dt, ...
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
    cost = max(abs(x(1)),abs(x(2)));
end

function cost = terminalcosts(t, x)
    cost = 0.0;
end

function [c,ceq] = constraints(t, x, u)
    global bound;
    c(1) = x(1)-bound;
    ceq = [];
end

function [c,ceq] = terminalconstraints(t, x)
    global bound;
    c(1) = x(1)-bound;
    ceq = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = -1;
    ub  =  1;
end

function y = system_dt(t, x, u, T)
    y(1) = ( -(x(1)^2+x(2)^2)*u(1)+x(1) )/...
           ( 1+(x(1)^2+x(2)^2)*u(1)^2-2*x(1)*u(1) );
    y(2) = x(2)/...
           ( 1+(x(1)^2+x(2)^2)*u(1)^2-2*x(1)*u(1) );
end

function dx = system_ct(t, x, u, T)
    dx = zeros(2,1);
    dx(1) = ( x(1)^2-x(2)^2 )*u(1);
    dx(2) = 2*x(1)*x(2)*u(1);
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
    global bound;
    [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
    figure(1);
        title('x_1/x_2 closed loop trajectory');
        xlabel('x_1');
        ylabel('x_2');
        grid on;
        hold on;
        for k=0.5:0.1:1.2
            plot(0.5*k*sin(0:0.01:2*pi),0.5*k*cos(0:0.01:2*pi)+0.5*k,'k');
            plot(0.5*k*sin(0:0.01:2*pi),0.5*k*cos(0:0.01:2*pi)-0.5*k,'k');
        end
        plot([bound bound], [-1.2 1.2], 'k', 'LineWidth', 2);
        plot(x_intermediate(:,1),x_intermediate(:,2),'-or', ...
             'MarkerFaceColor','r');
        axis([-1.2 1.2 -1.2 1.2]);
        axis square;

    figure(2);
        title(['x_1 and x_2 closed loop trajectory']);
        xlabel('n');
        ylabel('x_1(n), x_2(n)');
        grid on;
        hold on;
        plot(t_intermediate,x_intermediate(:,1),'-ok');
        plot(t_intermediate,x_intermediate(:,2),'-ok');
        axis([0 30 -0.5 1.1]);
        axis square;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
