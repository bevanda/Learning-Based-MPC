clear;
BEGIN_ACADO; %start ACADO
    acadoSet('problemname', 'mgcm');
    acadoSet('results_to_file',false);       % don't write results to file
    %% Initial setting of the model
    DifferentialState x1 x2; 
    Control u;
    
    B = 1;
    gamma = 1;
    x10 = 0;

    f = acado.DifferentialEquation();
    
    f.add(dot(x1) == 1/B*(x2-gamma*((x1+x10)^(1/2)-(x10)^(1/2))));
    f.add(dot(x2) == B*(-1*x2^3-1*x2^2-1*x2-x1-u));

    %% Initial setting of the optimisation problem
    N = 15; %prediction horizon
    sampling_time = 0.2;
    
    Q = 1e2*eye(3); 
    Q(1,1) = 1e2; 
    h = {u x1 x2}; 
    r = [0,0,0];
    
    ocp = acado.OCP(0.0, N*sampling_time, 20) ;
    ocp.minimizeLSQ(Q , h , r) ;
    ocp.subjectTo(f);
    ocp.subjectTo(0 <= x1 <= 1); 
    ocp.subjectTo(1.1875 <= x2 <= 2.1547);
    ocp.subjectTo(0.1547 <= u <= 2.1547);
    

       %% Setting up the process
    identity = acado.OutputFcn();
    
    dynamicSystem = acado.DynamicSystem(f, identity);
    
    process = acado.Process(dynamicSystem, 'INT_RK45');
    
    %% Setting up MPC controller
    Duration = 10;
    
    algo = acado.RealTimeAlgorithm(ocp, sampling_time);
    
    algo.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
    algo.set( 'MAX_NUM_ITERATIONS',          5                 );
    algo.set( 'INFEASIBLE_QP_HANDLING',      'YES'             );
    
%     algo.set('KKT_TOLERANCE', 1e-5);
%     algo.set('INTEGRATOR_TYPE', 'INT_RK45') 
%     algo.set('HESSIAN_APPROXIMATION', 'GAUSS_NEWTON');
%     algo.set('DISCRETIZATION_TYPE', 'COLLOCATION');

    reference = acado.StaticReferenceTrajectory();
    
    
    controller = acado.Controller(algo, reference);
    
    
    sim = acado.SimulationEnvironment(0.0, Duration, process, controller );
    
    
    r = [0.5,0.3,0.0];
    
    sim.init(r); 
    
    END_ACADO;       % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.
%% Run
out = mgcm_RUN();