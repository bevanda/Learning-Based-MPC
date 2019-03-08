function example_8_40iii

% Example 8.40(iii) from the second edition of L. Gruene, J. Pannek, 
% Nonlinear Model Predictive Control. Theory and Algorithms. Springer, 2017

numfig = 1;
for N=5:5:10

    [tc,xc,uc] = mpc_example_economic_tc(1.9,N,1,2,30,0);

    [tu,xu,uu] = mpc_example_economic(1.9,N,1,2,30,0);

    Kmin=1;
    Kmax=25;
    Kstep=1;

    for K = Kmin:Kstep:Kmax
        JKc(K) = sum(uc(1:K).^2);
        JKu(K) = sum(uu(1:K).^2);
    end

    figure(numfig);
    axes('FontSize',16)
    plot(Kmin:Kstep:Kmax,JKu,'kx--',Kmin:Kstep:Kmax,JKc,'ko-','LineWidth',1)
    hlegend = legend('no terminal conditions', 'terminal conditions','Location','SouthEast');
    set(hlegend,'Interpreter','latex','FontSize',16)
    ylb = sprintf('$$J_K^{cl}(1.9, \\mu_{%d})$$',N);
    latexxlabel('$$K$$',18)
    latexylabel(ylb,18)
    numfig = numfig + 1;
end;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main routine setting up the NMPC Problem with terminal conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [tt,xx,uu] = mpc_example_economic_tc(iv,N,T,constr,maxiter,myprint)

    addpath('./nmpcroutine');

    mpciterations = maxiter;
    t0            = 0.0;
    x0            = [iv];
    u0            = zeros(1,N);
    tol_opt       = 1e-12;
    opt_option    = 0;
    iprint        = 5;
    type          = 'difference equation';
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;

    global gconstr;
    gconstr = constr;
    
    global gN;
    gN = N;
    
    global fignum;
    if (myprint >=2)
        fignum = myprint-1;
        myprint = 2;
    end;
    
    global gprint;
    gprint = myprint;
    
    global guu;
    guu = [];
    
    global gtt;
    gtt = [ t0 ];
    
    global gxx;
    gxx = [ x0 ];
    
    nmpc(@runningcosts_tc, @terminalcosts_tc, @constraints_tc, ...
         @terminalconstraints_tc, @linearconstraints_tc, @system_tc, ...
         mpciterations, N, T, t0, x0, u0, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, ...
         rtol_ode_sim, iprint, @printHeader_tc, ...
         @printClosedloopData_tc, @plotTrajectories_tc);

    tt = gtt;
    xx = gxx;
    uu = guu;
    
    rmpath('./nmpcroutine');
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts_tc(t, x, u)
    cost = u(1)^2;
end

function cost = terminalcosts_tc(t, x)
    cost = 0.0;
end

function [c,ceq] = constraints_tc(t, x, u)
    global gconstr;
    
    c(1) = x(1) - gconstr;
    c(2) = - x(1) - gconstr;
    ceq = [];
end

function [c,ceq] = terminalconstraints_tc(t, x)
    global gconstr;
    
    c(1) = x(1) - gconstr;
    c(2) = - x(1) - gconstr;
    ceq = [x(1)];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints_tc(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb = -3.0;
    ub = 3.0;
end

function y = system_tc(t, x, u, T)
    y(1) = 2.0.*x(1)+T*u(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of output format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printHeader_tc()
    fprintf('   k  |      u(k)        x(1)        x(2)     Time\n');
    fprintf('--------------------------------------------------\n');
end

function printClosedloopData_tc(mpciter, u, x, t_Elapsed)
    
    fprintf(' %3d  | %+11.6f %+11.6f  %+6.3f', mpciter, ...
            u(1,1), x(1), t_Elapsed);
       
end

function plotTrajectories_tc(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    global gN;
    global fignum;
    global gprint;
    global gtt;
    global gxx;
    global guu;
    global gconstr;
    
    global count;
    
    gtt = [gtt; t0+T];
    
    if (t0==0)
        count = 0;
    end;
    
    if (gprint>=1)
        figure(fignum);
        axis([0, 25, -0.2, 2.2])
        hold on;  
        xlabel('n');
        ylabel('x(n)');
        grid on;
    end;
    x = x0;
    for k=1:gN
        [x,t_intermediate, x_intermediate] = ...
            dynamic(system, T, t0+k-1, x, u(:,k), ...
                         atol_ode, rtol_ode, type);
        if (k==1)
            if (gprint>=1)
                clt = t_intermediate;
                clx = x_intermediate;
            end;
            gxx = [gxx; x];
            guu = [guu, u(1,1)];
        end
        if (gprint>=2)
            plot(t_intermediate,x_intermediate(:,1),'--k',...
                 'linewidth',2);
        end 
    end
     if (gprint>=1)
        plot(clt,clx,'-k','linewidth',2);
    end;
    
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main routine setting up the NMPC Problem without terminal conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [tt,xx,uu,val] = mpc_example_economic(iv,N,T,constr,mpciter,myprint)

    addpath('./nmpcroutine');

    mpciterations = mpciter;
    t0            = 0.0;
    x0            = iv;
    u0            = zeros(1,N);
    tol_opt       = 1e-12;
    opt_option    = 0;
    iprint        = 5;
    type          = 'difference equation';
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;

    global gconstr;
    gconstr = constr;
    
    global gN;
    gN = N;
    
    global fignum;
    if (myprint >=2)
        fignum = myprint-1;
        myprint = 2;
    end;
    
    global gprint;
    gprint = myprint;
    
    global guu;
    guu = [];
    
    global gtt;
    gtt = [ t0 ];
    
    global gxx;
    gxx = [ x0 ];
    
    global gval;
    
    
    
    nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system, ...
         mpciterations, N, T, t0, x0, u0, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint, @printHeader, @printClosedloopData, @plotTrajectories);

    tt = gtt;
    xx = gxx;
    uu = guu;
    val = gval;
    
    rmpath('./nmpcroutine');
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)
    cost = u(1)^2;
end

function cost = terminalcosts(t, x)
    cost = 0.0;
end

function [c,ceq] = constraints(t, x, u)
    global gconstr;
    
    c(1) = x(1) - gconstr;
    c(2) = - x(1) - gconstr;
    ceq = [];
end

function [c,ceq] = terminalconstraints(t, x)
    global gconstr;
    
    c(1) = x(1) - gconstr;
    c(2) = - x(1) - gconstr;
    ceq = [];
    %ceq = [x(1)-0.5];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb = -3.0;
    ub = 3.0;
end

function y = system(t, x, u, T)
    y(1) = 2.0.*x(1)+T*u(1);
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
    
    fprintf(' %3d  | %+11.6f %+11.6f  %+6.3f', mpciter, ...
            u(1,1), x(1), t_Elapsed);
       
end

function plotTrajectories(dynamic, system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    global gN;
    global gprint;
    global fignum;
    global gtt;
    global gxx;
    global guu;
    global gval;
    
    gtt = [gtt; t0+T];
    
    if (gprint>=1)
        figure(fignum);
        hold on;  
    end;
        x = x0;
        gval = 0;
        for k=1:gN
            gval = gval + runningcosts(t0+k-1,x,u(:,k));
            [x,t_intermediate, x_intermediate] = ...
                dynamic(system, T, t0+k-1, x, u(:,k), ...
                             atol_ode, rtol_ode, type);
            if (k==1)
                if (gprint>=1)
                    plot(t_intermediate,x_intermediate(:,1),'-k');
                end;
                gxx = [gxx; x];
                guu = [guu, u(1,1)];
            else
                if (gprint>=2)
                    plot(t_intermediate,x_intermediate(:,1),'--k');
                end;
            end 
        end
        gval = gval/gN;
        if (gprint>=1)
            if (gprint>=2)
                title(['open loop trajectories (dashed) and closed loop trajectory (solid)']);
            else
                title(['closed loop trajectory (solid)']);
            end;
            xlabel('n');
            ylabel('x(n)');
            grid on;
        end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions for fancy captions and labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  s=latextitle(string,size)

set(get(gca,'title'),...
    'Interpreter','latex',...
	'String',string,...
    'FontName','times',...
	'FontSize',size);
end


function  s=latexxlabel(string,size)

set(get(gca,'Xlabel'),...
    'Interpreter','latex',...
	'String',string,...
    'FontName','times',...
	'FontSize',size);

end

function  s=latexylabel(string,size)

set(get(gca,'Ylabel'),...
    'Interpreter','latex',...
	'String',string,...
    'FontName','times',...
	'FontSize',size);
end


