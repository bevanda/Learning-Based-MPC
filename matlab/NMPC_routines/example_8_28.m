function example_8_28

% Example 8.28 from the second edition of L. Gruene, J. Pannek, Nonlinear 
% Model Predictive Control. Theory and Algorithms. Springer, 2017

Nmin = 2;
Nmax = 15;
iter = 30;
clear d;

for N = Nmin:Nmax
    [t x u] = mpc_example_economic(1.9,N,1,2,iter,0);
    d(N-Nmin+1) = u(iter)^2;
    fprintf('d(%d)=%f\n',N,d(N-Nmin+1));
end

figure;
semilogy(Nmin:Nmax,d,'--*k','linewidth',2);
latexxlabel('$$N$$',18);
latexylabel('$$\overline J_{\infty}^{cl}(1.9,\mu_N)$$',18);
axis([1.5,15.5,1e-8,1e1]);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main routine setting up the NMPC Problem without terminal conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


