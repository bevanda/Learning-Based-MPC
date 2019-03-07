%Discrete Optimal Control Problem
 clear;

BEGIN_ACADO;        %start ACADO

    acadoSet('problemname','swing_up');      % set problem name
    acadoSet('results_to_file',false);       % don't write results to file
    

    DifferentialState z;                     % cart position
    DifferentialState z_dot;                 % cart velocity 
    DifferentialState theta;                 % angle  
    DifferentialState theta_dot;             % angular velocity
    
    DifferentialState F;                     % control force
    
    
    Control F_dot;                           % control change
    
    
    mCart = 1;                               % cart mass
    mPend = 1;                               % pendulum mass
    g = 9.81;                                % gravity of earth
    L = 0.5;                                 % pendulum length
    Kd = 10;                                 % cart damping

    
    
    
    %% Differential Equations
    f = acado.DifferentialEquation();
    

    f.add( dot(z) == z_dot );
    f.add( dot(z_dot) == (F - Kd*z_dot - mPend*L*theta_dot^2*sin(theta) + mPend*g*sin(theta)*cos(theta)) / (mCart + mPend*sin(theta)^2) );
    f.add( dot(theta) == theta_dot );
    f.add( dot(theta_dot) == ((F - Kd*z_dot - mPend*L*theta_dot^2*sin(theta))*cos(theta)/(mCart + mPend) + g*sin(theta)) / (L - mPend*L*cos(theta)^2/(mCart + mPend)) );
    f.add( dot(F) == F_dot );


    
    %% Optimal Control Problem 
    N = 25;
    Sampling_time = 0.1;
    
    ocp = acado.OCP(0.0,N*Sampling_time, 20);
    
    h = { z, z_dot, theta, theta_dot, F_dot };
    
    Q = diag([100,1,100,1,0.01]);
    
    r = [0,0,0,0,0];
    
    
    ocp.minimizeLSQ( Q,h,r ); 
    
    
    % constraints
    ocp.subjectTo( f );
    
    ocp.subjectTo( -10 <= z <= 10 );
    ocp.subjectTo( -100 <= F <= 100 );
    
    
    %% Setting up the process
    identity = acado.OutputFcn();
    
    dynamicSystem = acado.DynamicSystem(f, identity);
    
    process = acado.Process(dynamicSystem, 'INT_RK45');
    
    %% Setting up MPC controller
    Duration = 10;
    
    algo = acado.RealTimeAlgorithm(ocp, Sampling_time);
    
    algo.set( 'DISCRETIZATION_TYPE',         'COLLOCATION' );
    algo.set( 'MAX_NUM_ITERATIONS',          5                );
    algo.set( 'INFEASIBLE_QP_HANDLING',      'YES'             );
    
    
    reference = acado.StaticReferenceTrajectory();
    
    
    controller = acado.Controller( algo, reference );
    
    
    sim = acado.SimulationEnvironment (0.0, Duration, process, controller );
    
    
    r = [0,0,-pi,0,0];
    
    sim.init(r);
    
    
    END_ACADO;           % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.
    
