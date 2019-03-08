function [t, x, u] = nmpc(runningcosts, terminalcosts, ...
              constraints, terminalconstraints, ...
              linearconstraints, system, ...
              mpciterations, N, T, tmeasure, xmeasure, u0, ...
              varargin)
% nmpc(runningcosts, terminalcosts, constraints, ...
%      terminalconstraints, linearconstraints, system, ...
%      mpciterations, N, T, tmeasure, xmeasure, u0, ...
%      tol_opt, opt_option, ...
%      type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
%      iprint, printHeader, printClosedloopData, plotTrajectories)
% Computes the closed loop solution for the NMPC problem defined by
% the functions
%   runningcosts:         evaluates the running costs for state and control
%                         at one sampling instant.
%                         The function returns the running costs for one
%                         sampling instant.
%          Usage: [cost] = runningcosts(t, x, u)
%                 with time t, state x and control u
%   terminalcosts:        evaluates the terminal costs for state at the end
%                         of the open loop horizon.
%                         The function returns value of the terminal costs.
%          Usage: cost = terminalcosts(t, x)
%                 with time t and state x
%   constraints:          computes the value of the restrictions for a
%                         sampling instance provided the data t, x and u
%                         given by the optimization method.
%                         The function returns the value of the
%                         restrictions for a sampling instance separated
%                         for inequality restrictions c and equality
%                         restrictions ceq.
%          Usage: [c,ceq] = constraints(t, x, u)
%                 with time t, state x and control u
%   terminalconstraints:  computes the value of the terminal restrictions
%                         provided the data t, x and u given by the
%                         optimization method.
%                         The function returns the value of the
%                         terminal restriction for inequality restrictions
%                         c and equality restrictions ceq.
%          Usage: [c,ceq] = terminalconstraints(t, x)
%                 with time t and state x
%   linearconstraints:    sets the linear constraints of the discretized
%                         optimal control problem. This is particularly
%                         useful to set control and state bounds.
%                         The function returns the required matrices for
%                         the linear inequality and equality constraints A
%                         and Aeq, the corresponding right hand sides b and
%                         beq as well as the lower and upper bound of the
%                         control.
%          Usage: [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
%                 with time t, state x and control u
%   system:               evaluates the difference equation describing the
%                         process given time t, state vector x and control
%                         u.
%                         The function returns the state vector x at the
%                         next time instant.
%          Usage: [y] = system(t, x, u, T)
%                 with time t, state x, control u and sampling interval T
% for a given number of NMPC iteration steps (mpciterations). For
% the open loop problem, the horizon is defined by the number of
% time instances N and the sampling time T. Note that the dynamic
% can also be the solution of a differential equation. Moreover, the
% initial time tmeasure, the state measurement xmeasure and a guess of
% the optimal control u0 are required.
%
% Arguments:
%   mpciterations:  Number of MPC iterations to be performed
%   N:              Length of optimization horizon
%   T:              Sampling interval
%   tmeasure:       Time measurement of initial value
%   xmeasure:       State measurement of initial value
%   u0:             Initial guess of open loop control
%
% Optional arguments:
%       iprint:     = 0  Print closed loop data(default)
%                   = 1  Print closed loop data and errors of the
%                        optimization method
%                   = 2  Print closed loop data and errors and warnings of
%                        the optimization method
%                   >= 5 Print closed loop data and errors and warnings of
%                        the optimization method as well as graphical
%                        output of closed loop state trajectories
%                   >=10 Print closed loop data and errors and warnings of
%                        the optimization method with error and warning
%                        description
%   printHeader:         Clarifying header for selective output of closed
%                        loop data, cf. printClosedloopData
%   printClosedloopData: Selective output of closed loop data
%   plotTrajectories:    Graphical output of the trajectories, requires
%                        iprint >= 4
%       tol_opt:         Tolerance of the optimization method
%       opt_option: = 0: Active-set method used for optimization (default)
%                   = 1: Interior-point method used for optimization
%                   = 2: Trust-region reflective method used for
%                        optimization
%   type:                Type of dynamic, either difference equation or
%                        differential equation can be used
%    atol_ode_real:      Absolute tolerance of the ODE solver for the
%                        simulated process
%    rtol_ode_real:      Relative tolerance of the ODE solver for the
%                        simulated process
%    atol_ode_sim:       Absolute tolerance of the ODE solver for the
%                        simulated NMPC prediction
%    rtol_ode_sim:       Relative tolerance of the ODE solver for the
%                        simulated NMPC prediction
%
% Internal Functions:
%   measureInitialValue:          measures the new initial values for t0
%                                 and x0 by adopting values computed by
%                                 method applyControl.
%                                 The function returns new initial state
%                                 vector x0 at sampling instant t0.
%   applyControl:                 applies the first control element of u to
%                                 the simulated process for one sampling
%                                 interval T.
%                                 The function returns closed loop state
%                                 vector xapplied at sampling instant
%                                 tapplied.
%   shiftHorizon:                 applies the shift method to the open loop
%                                 control in order to ease the restart.
%                                 The function returns a new initial guess
%                                 u0 of the control.
%   solveOptimalControlProblem:   solves the optimal control problem of the
%                                 horizon N with sampling length T for the
%                                 given initial values t0 and x0 and the
%                                 initial guess u0 using the specified
%                                 algorithm.
%                                 The function returns the computed optimal
%                                 control u, the corresponding value of the
%                                 cost function V as well as possible exit
%                                 flags and additional output of the
%                                 optimization method.
%   costfunction:                 evaluates the cost function of the
%                                 optimal control problem over the horizon
%                                 N with sampling time T for the current
%                                 data of the optimization method t0, x0
%                                 and u.
%                                 The function return the computed cost
%                                 function value.
%   nonlinearconstraints:         computes the value of the restrictions
%                                 for all sampling instances provided the
%                                 data t0, x0 and u given by the
%                                 optimization method.
%                                 The function returns the value of the
%                                 restrictions for all sampling instances
%                                 separated for inequality restrictions c
%                                 and equality restrictions ceq.
%   computeOpenloopSolution:      computes the open loop solution over the
%                                 horizon N with sampling time T for the
%                                 initial values t0 and x0 as well as the
%                                 control u.
%                                 The function returns the complete open
%                                 loop solution over the requested horizon.
%   dynamic:                      evaluates the dynamic of the system for
%                                 given initial values t0 and x0 over the
%                                 interval [t0, tf] using the control u.
%                                 The function returns the state vector x
%                                 at time instant tf as well as an output
%                                 of all intermediate evaluated time
%                                 instances.
%   printSolution:                prints out information on the current MPC
%                                 step, in particular state and control
%                                 information as well as required computing
%                                 times and exitflags/outputs of the used
%                                 optimization method. The flow of
%                                 information can be controlled by the
%                                 variable iprint and the functions
%                                 printHeader, printClosedloopData and
%                                 plotTrajectories.
%
% Version of May 30, 2011, in which a bug appearing in the case of 
% multiple constraints has been fixed
%
% (C) Lars Gruene, Juergen Pannek 2011

    if (nargin>=13)
        tol_opt = varargin{1};
    else
        tol_opt = 1e-6;
    end;
    if (nargin>=14)
        opt_option = varargin{2};
    else
        opt_option = 0;
    end;
    if (nargin>=15)
        if ( strcmp(varargin{3}, 'difference equation') || ...
                strcmp(varargin{3}, 'differential equation') )
            type = varargin{3};
        else
            fprintf([' Wrong input for type of dynamic: use either ', ...
                '"difference equation" or "differential equation".']);
        end
    else
        type = 'difference equation';
    end;
    if (nargin>=16)
        atol_ode_real = varargin{4};
    else
        atol_ode_real = 1e-8;
    end;
    if (nargin>=17)
        rtol_ode_real = varargin{5};
    else
        rtol_ode_real = 1e-8;
    end;
    if (nargin>=18)
        atol_ode_sim = varargin{6};
    else
        atol_ode_sim = atol_ode_real;
    end;
    if (nargin>=19)
        rtol_ode_sim = varargin{7};
    else
        rtol_ode_sim = rtol_ode_real;
    end;
    if (nargin>=20)
        iprint = varargin{8};
    else
        iprint = 0;
    end;
    if (nargin>=21)
        printHeader = varargin{9};
    else
        printHeader = @printHeaderDummy;
    end;
    if (nargin>=22)
        printClosedloopData = varargin{10};
    else
        printClosedloopData = @printClosedloopDataDummy;
    end;
    if (nargin>=23)
        plotTrajectories = varargin{11};
    else
        plotTrajectories = @plotTrajectoriesDummy;
    end;

    % Determine MATLAB Version and
    % specify and configure optimization method
    vs = version('-release');
    vyear = str2num(vs(1:4));
    if (vyear <= 2007)
        fprintf('MATLAB version R2007 or earlier detected\n');
        if ( opt_option == 0 )
            options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 2000,...
                'LargeScale', 'off',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1);
        elseif ( opt_option == 1 )
            error('nmpc:WrongArgument', '%s\n%s', ...
                  'Interior point method not supported in MATLAB R2007', ...
                  'Please use opt_option = 0 or opt_option = 2');
        elseif ( opt_option == 2 )
             options = optimset('Display','off',...
                 'TolFun', tol_opt,...
                 'MaxIter', 2000,...
                 'LargeScale', 'on',...
                 'Hessian', 'off',...
                 'MaxPCGIter', max(1,floor(size(u0,1)*size(u0,2)/2)),...
                 'PrecondBandWidth', 0,...
                 'TolPCG', 1e-1);
        end
    else
        fprintf('MATLAB version R2008 or newer detected\n');
        if ( opt_option == 0 )
            options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 10000,...
                'Algorithm', 'active-set',...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolConSQP', 1e-6);
        elseif ( opt_option == 1 )
            options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 2000,...
                'Algorithm', 'interior-point',...
                'AlwaysHonorConstraints', 'bounds',...
                'FinDiffType', 'forward',...
                'HessFcn', [],...
                'Hessian', 'bfgs',...
                'HessMult', [],...
                'InitBarrierParam', 0.1,...
                'InitTrustRegionRadius', sqrt(size(u0,1)*size(u0,2)),...
                'MaxProjCGIter', 2*size(u0,1)*size(u0,2),...
                'ObjectiveLimit', -1e20,...
                'ScaleProblem', 'obj-and-constr',...
                'SubproblemAlgorithm', 'cg',...
                'TolProjCG', 1e-2,...
                'TolProjCGAbs', 1e-10);
        %                       'UseParallel','always',...
        elseif ( opt_option == 2 )
            options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 2000,...
                'Algorithm', 'trust-region-reflective',...
                'Hessian', 'off',...
                'MaxPCGIter', max(1,floor(size(u0,1)*size(u0,2)/2)),...
                'PrecondBandWidth', 0,...
                'TolPCG', 1e-1);
        end
    end

    warning off all
    t = [];
    x = [];
    u = [];
    % Start of the NMPC iteration
    mpciter = 0;
    while(mpciter < mpciterations)
        % Step (1) of the NMPC algorithm:
        %   Obtain new initial value
        [t0, x0] = measureInitialValue ( tmeasure, xmeasure );
        % Step (2) of the NMPC algorithm:
        %   Solve the optimal control problem
        t_Start = tic;
        [u_new, V_current, exitflag, output] = solveOptimalControlProblem ...
            (runningcosts, terminalcosts, constraints, ...
            terminalconstraints, linearconstraints, system, ...
            N, t0, x0, u0, T, ...
            atol_ode_sim, rtol_ode_sim, tol_opt, options, type);
        t_Elapsed = toc( t_Start );
        %   Print solution
        if ( iprint >= 1 )
            printSolution(system, printHeader, printClosedloopData, ...
                          plotTrajectories, mpciter, T, t0, x0, u_new, ...
                          atol_ode_sim, rtol_ode_sim, type, iprint, ...
                          exitflag, output, t_Elapsed);
        end
        %   Store closed loop data
        t = [ t; tmeasure ];
        x = [ x; xmeasure ];
        u = [ u; u_new(:,1) ];
        %   Prepare restart
        u0 = shiftHorizon(u_new);
        % Step (3) of the NMPC algorithm:
        %   Apply control to process
        [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_new, ...
            atol_ode_real, rtol_ode_real, type);
        mpciter = mpciter+1;
    end
end

function [t0, x0] = measureInitialValue ( tmeasure, xmeasure )
    t0 = tmeasure;
    x0 = xmeasure;
end

function [tapplied, xapplied] = applyControl(system, T, t0, x0, u, ...
                                atol_ode_real, rtol_ode_real, type)
    xapplied = dynamic(system, T, t0, x0, u(:,1), ...
                       atol_ode_real, rtol_ode_real, type);
    tapplied = t0+T;
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function [u, V, exitflag, output] = solveOptimalControlProblem ...
    (runningcosts, terminalcosts, constraints, terminalconstraints, ...
    linearconstraints, system, N, t0, x0, u0, T, ...
    atol_ode_sim, rtol_ode_sim, tol_opt, options, type)
    x = zeros(N+1, length(x0));
    x = computeOpenloopSolution(system, N, T, t0, x0, u0, ...
                                atol_ode_sim, rtol_ode_sim, type);

    % Set control and linear bounds
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    for k=1:N
        [Anew, bnew, Aeqnew, beqnew, lbnew, ubnew] = ...
               linearconstraints(t0+k*T,x(k,:),u0(:,k));
        A = blkdiag(A,Anew);
        b = [b, bnew];
        Aeq = blkdiag(Aeq,Aeqnew);
        beq = [beq, beqnew];
        lb = [lb, lbnew];
        ub = [ub, ubnew];
    end

    % Solve optimization problem
    [u, V, exitflag, output] = fmincon(@(u) costfunction(runningcosts, ...
        terminalcosts, system, N, T, t0, x0, ...
        u, atol_ode_sim, rtol_ode_sim, type), u0, A, b, Aeq, beq, lb, ...
        ub, @(u) nonlinearconstraints(constraints, terminalconstraints, ...
        system, N, T, t0, x0, u, ...
        atol_ode_sim, rtol_ode_sim, type), options);
end

function cost = costfunction(runningcosts, terminalcosts, system, ...
                    N, T, t0, x0, u, ...
                    atol_ode_sim, rtol_ode_sim, type)
    cost = 0;
    x = zeros(N+1, length(x0));
    x = computeOpenloopSolution(system, N, T, t0, x0, u, ...
                                atol_ode_sim, rtol_ode_sim, type);
    for k=1:N
        cost = cost+runningcosts(t0+k*T, x(k,:), u(:,k));
    end
    cost = cost+terminalcosts(t0+(N+1)*T, x(N+1,:));
end

function [c,ceq] = nonlinearconstraints(constraints, ...
    terminalconstraints, system, ...
    N, T, t0, x0, u, atol_ode_sim, rtol_ode_sim, type)
    x = zeros(N+1, length(x0));
    x = computeOpenloopSolution(system, N, T, t0, x0, u, ...
                                atol_ode_sim, rtol_ode_sim, type);
    c = [];
    ceq = [];
    for k=1:N
        [cnew, ceqnew] = constraints(t0+k*T,x(k,:),u(:,k));
        c = [c cnew];
        ceq = [ceq ceqnew];
    end
    [cnew, ceqnew] = terminalconstraints(t0+(N+1)*T,x(N+1,:));
    c = [c cnew];
    ceq = [ceq ceqnew];
end

function x = computeOpenloopSolution(system, N, T, t0, x0, u, ...
                                     atol_ode_sim, rtol_ode_sim, type)
    x(1,:) = x0;
    for k=1:N
        x(k+1,:) = dynamic(system, T, t0, x(k,:), u(:,k), ...
                             atol_ode_sim, rtol_ode_sim, type);
    end
end

function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
             x0, u, atol_ode, rtol_ode, type)
    if ( strcmp(type, 'difference equation') )
        x = system(t0, x0, u, T);
        x_intermediate = [x0; x];
        t_intermediate = [t0, t0+T];
    elseif ( strcmp(type, 'differential equation') )
        options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
        [t_intermediate,x_intermediate] = ode45(system, ...
            [t0, t0+T], x0, options, u);
        x = x_intermediate(size(x_intermediate,1),:);
    end
end

function printSolution(system, printHeader, printClosedloopData, ...
             plotTrajectories, mpciter, T, t0, x0, u, ...
             atol_ode, rtol_ode, type, iprint, exitflag, output, t_Elapsed)
    if (mpciter == 0)
        printHeader();
    end
    printClosedloopData(mpciter, u, x0, t_Elapsed);
    switch exitflag
        case -2
        if ( iprint >= 1 && iprint < 10 )
            fprintf(' Error F\n');
        elseif ( iprint >= 10 )
            fprintf(' Error: No feasible point was found\n')
        end
        case -1
        if ( iprint >= 1 && iprint < 10 )
            fprintf(' Error OT\n');
        elseif ( iprint >= 10 )
            fprintf([' Error: The output function terminated the',...
                     ' algorithm\n'])
        end
        case 0
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning IT\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Number of iterations exceeded',...
                     ' options.MaxIter or number of function',...
                     ' evaluations exceeded options.FunEvals\n'])
        end
        case 1
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' \n');
        elseif ( iprint >= 10 )
            fprintf([' First-order optimality measure was less',...
                     ' than options.TolFun, and maximum constraint',...
                     ' violation was less than options.TolCon\n'])
        end
        case 2
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning TX\n');
        elseif ( iprint >= 10 )
            fprintf(' Warning: Change in x was less than options.TolX\n')
        end
        case 3
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning TJ\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Change in the objective function',...
                     ' value was less than options.TolFun\n'])
        end
        case 4
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning S\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Magnitude of the search direction',...
                     ' was less than 2*options.TolX and constraint',...
                     ' violation was less than options.TolCon\n'])
        end
        case 5
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning D\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Magnitude of directional derivative',...
                     ' in search direction was less than',...
                     ' 2*options.TolFun and maximum constraint',...
                     ' violation was less than options.TolCon\n'])
        end
    end
    if ( iprint >= 5 )
        plotTrajectories(@dynamic, system, T, t0, x0, u, atol_ode, rtol_ode, type)
    end
end

function printHeaderDummy(varargin)
end

function printClosedloopDataDummy(varargin)
end

function plotTrajectoriesDummy(varargin)
end







