function mgcm_nmpc

    addpath('./nmpcroutine');
    clearvars;
    close all;

    mpciterations = 100;
    N             = 3;
    T             = 0.01;
    tmeasure      = 0.0;
    xmeasure      = [0.150000000000000,1.28750000000000,1.15470000000000,0]; % init state measurement
    u0            = 0.2*ones(1,N);
    tol_opt       = 1e-8;
    opt_option    = 0;
    iprint        = 10;
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
    Q = eye(4);  R = eye(1);
    %cost = x(1)^2+x(2)^2+x(3)^2+x(4)^2 +u(1)^2;
    cost = x*Q*x'+u*R*u';
     
end

function cost = terminalcosts(t, x)
    % P solution to a lyapunov equation
    % T punishment for deviation to the steady state
    T = 1e3;
    P = [8292.85866674464,3829.78519283510,-1719.52687620407,-16.5895240771607;...
        3829.78519283510,1830.09184563230,-828.404573677948,-7.99393347184673;...
        -1719.52687620407,-828.404573677948,493.014176560514,4.79962481767885; ...
        -16.5895240771607,-7.99393347184673,4.79962481767885,1.04998423860038];
%     cost = x*P*x';
    cost = 0.0;
end

function [c,ceq] = constraints(t, x, u)
    % Constraints in the form --->  c() <= 0, ceq()== 0
    mflow_min=0; mflow_max=1;
    prise_min=1.1875; prise_max=2.1875;
    throttle_min=0.1547; throttle_max=2.1547;
    throttle_rate_min=-20; throttle_rate_max=20;
    u_min=0.1547;u_max=2.1547;
    
    c(1) =  x(1)-mflow_max; 
    c(2) = -x(1)+mflow_min;
    c(3) =  x(2)-prise_max;
    c(4) = -x(2)+prise_min;
    c(5) = -x(3)-throttle_max;
    c(6) =  x(3)+throttle_min;
    c(7) =  x(4)-throttle_rate_max;
    c(8) = -x(4)+throttle_rate_min;
    c(9) = u(1)-u_max;
    c(10) = -u(1)+u_min;
    ceq  = [];
end

function [c,ceq] = terminalconstraints(t, x, u)
    % Constraints in the form --->  c() <= 0, ceq()== 0
    mflow_min=0; mflow_max=1;
    prise_min=1.1875; prise_max=2.1875;
    throttle_min=0.1547; throttle_max=2.1547;
    throttle_rate_min=-20; throttle_rate_max=20;
    u_min=0.1547; u_max=2.1547;
    
    c(1) =  x(1)-mflow_max; 
    c(2) = -x(1)+mflow_min;
    c(3) =  x(2)-prise_max;
    c(4) = -x(2)+prise_min;
    c(5) = -x(3)-throttle_max;
    c(6) =  x(3)+throttle_min;
    c(7) =  x(4)-throttle_rate_max;
    c(8) = -x(4)+throttle_rate_min;
    ceq  = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = [];
    ub  = [];
end

function y = system_dt(t, x, u, T)
% Discrete time state-space model of the Moore-Greitzer compressor model
    Ad = [1.01125000000000,0.0100000000000000,0,0;...
        0.0100000000000000,0.995555557627778,-0.0129903810567666,0;...
        0,0,1,0.0100000000000000;...
        0,0,-10,0.552786404500042];
    Bd = [0;0;0;10];
    K=[-3.0741 2.0957 0.1197 -0.0090];
    z = Ad*x'+Bd*u;
    y = z';
end

function dx = system_ct(t, x, u, T)
% Continous time state-space model of the Moore-Greitzer compressor model
    dx = zeros(4,1);
    wn = sqrt(1000); % resonant frequency
    zeta = 1/sqrt(2); % damping coefficient
    beta = 1; % constant >0
    x2_c = 0; % pressure constant
   
    
    dx(1) = x(2)+x2_c+1+3*(x(1)/2)-(x(1)^3/2); % mass flow rate
    dx(2) = (x(1)+1-x(3)*sqrt(x(2)))/(beta^2); % pressure rise rate
    dx(3) = x(4); % throttle opening rate
    dx(4) = -wn^2*x(3)-2*zeta*wn*x(4)+wn^2*u(1); % throttle opening acceleration

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of output format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printHeader()
    fprintf('   k  |      u(k)        x(1)        x(2)        x(3)      x(4)    Time\n');
    fprintf('----------------------------------------------------------------------\n');
end

function printClosedloopData(mpciter, u, x, t_Elapsed)
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+6.3f', mpciter, ...
            u(1,1), x(1), x(2), x(3), x(4), t_Elapsed);
end

function plotTrajectories(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
    figure(1);
        title('x_1/x_2 closed loop trajectory');
        xlabel('mass flow');
        ylabel('pressure rise');
        grid on;
        hold on;
%         plot(0, -0.5,'ok','MarkerSize',8);
        plot(x_intermediate(:,1),x_intermediate(:,2),'-ok');
%         axis([-1 1 -1 1]);
        axis square;

    figure(2);
        title(['x_1 and x_2 closed loop trajectory']);
        xlabel('n');
        ylabel('x_1(n), x_2(n)');
        grid on;
        hold on;
        size(t_intermediate)
        size(x_intermediate(:,1))
        plot(t_intermediate,x_intermediate(:,1),'-ok');
        plot(t_intermediate,x_intermediate(:,2),'-ok');
        axis([0 mpciterations -0.5 1]);
        axis square;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
