function example_8_1

% Example 8.1 from the second edition of L. Gruene, J. Pannek, Nonlinear
% Model Predictive Control. Theory and Algorithms. Springer, 2017

[t,x,u]=mpc_example_economic_tc(1.9,3,1,2,30,2);
figure(1);
latexxlabel('$$k$$',18);
latexylabel('$$x_{\mu_N}(k)$$ (solid) and $$x_{u^*_N}(\cdot)$$ (dashed)',18);
plot([0,35],[2,2],'k-');
plot([0,35],[0,0],'k--');
latextitle('N=3',18);
axis([0,25,-0.1,2.1]);


[t,x,u]=mpc_example_economic_tc(1.9,10,1,2,30,3);
figure(2);
latexxlabel('$$k$$',18);
latexylabel('$$x_{\mu_N}(k)$$ (solid) and $$x_{u^*_N}(\cdot)$$ (dashed)',18);
plot([0,40],[2,2],'k-');
plot([0,40],[0,0],'k--');
latextitle('N=10',18);
axis([0,40,-0.1,2.1]);

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
    
    nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system, ...
         mpciterations, N, T, t0, x0, u0, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint, @printHeader, @printClosedloopData, @plotTrajectories);

    tt = gtt;
    xx = gxx;
    uu = guu;
    
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
    ceq = [x(1)];
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

