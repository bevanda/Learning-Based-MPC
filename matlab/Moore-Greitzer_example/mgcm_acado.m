clear;
    BEGIN_ACADO; %start ACADO
    acadoSet('problemname', 'mgcm');
    acadoSet('results_to_file',false);       % don't write results to file
    %% Initial setting of the model
    DifferentialState x1 x2 x3 x4; 
    Control u;
    
    wn=sqrt(1000); % resonant frequency
    zeta=1/sqrt(2); % damping coefficient
    beta=1; % constant >0
    x2_c=0; % pressure constant
    
    mflow_min=0; mflow_max=1;
    prise_min=1.1875; prise_max=2.1875;
    throttle_min=0.1547; throttle_max=2.1547;
    throttle_rate_min=-20; throttle_rate_max=20;
    u_min=0.1547;u_max=2.1547;

    f = acado.DifferentialEquation();
    
    f.add(dot(x1) == -x2+x2_c+1+3*(x1/2)-(x1^3/2)); % mass flow rate
    f.add(dot(x2) == (x1+1-x3*sqrt(x2))/(beta^2)); % pressure rise rate
    f.add(dot(x3) == x4); % throttle opening rate
    f.add(dot(x4) == -wn^2*x3-2*zeta*wn*x4+wn^2*u) ;% throttle opening acceleration

    %% Initial setting of the optimisation problem
    duration = 600; %length of simulation
    N = 50; %prediction horizon
    sampling_time = 0.01;
    
    Q = eye(5);  
    h = {u x1 x2 x3 x4}; 
        %% init conditions without acado real time sim
    x10 = 0.5;
    x20 = 1.6875;
    x30 = 1.1547;
    x40 = 0;
    r = [x10,x20,x30,x40,x30];
    
    ocp = acado.OCP(0.0, duration) ;
    ocp.minimizeLSQ(Q, h, r);
    %% Constraints
    ocp.subjectTo(f);
%     ocp.subjectTo(mflow_min <= x1 <= mflow_max); 
%     ocp.subjectTo(prise_min <= x2 <= prise_max);
%     ocp.subjectTo(throttle_min <= x3 <= throttle_max); 
%     ocp.subjectTo(throttle_rate_min <= x4 <= throttle_rate_max)
%     ocp.subjectTo(u_min <= u <= u_max);
    

%     ocp.subjectTo('AT_START', x1==x10-0.35);
%     ocp.subjectTo('AT_START', x2==x20-0.4);
%     ocp.subjectTo('AT_START', x3==x30);
%     ocp.subjectTo('AT_START', x4==0);

    %% Setting up the process
%     identity = acado.OutputFcn();
%     
%     dynamicSystem = acado.DynamicSystem(f, identity);
%     
%     process = acado.Process(dynamicSystem, 'INT_RK45');
    
    %% Setting up MPC controller
    P = 5; %maximal number of SQP iterations 
%   algo = acado.RealTimeAlgorithm(ocp, sampling_time); % real-time algorithm
    algo = acado.OptimizationAlgorithm(ocp); 
    algo.set( 'DISCRETIZATION_TYPE', 'COLLOCATION' );
    algo.set( 'MAX_NUM_ITERATIONS', P);
%     algo.set( 'INFEASIBLE_QP_HANDLING', 'YES' );
    
    algo.set('KKT_TOLERANCE', 1e-6);
    algo.set('INTEGRATOR_TYPE', 'INT_RK78') 
    algo.set('HESSIAN_APPROXIMATION', 'GAUSS_NEWTON');

%     reference = acado.StaticReferenceTrajectory();
%     
%     
%     controller = acado.Controller(algo, reference);
%     
%     
%     sim = acado.SimulationEnvironment(0.0, Duration, process, controller );
%     
 
%     r = [x10-0.35, x20-0.4, x30, 0];  % init condition
%     sim.init(r); %set init values
    
    END_ACADO;       % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.
