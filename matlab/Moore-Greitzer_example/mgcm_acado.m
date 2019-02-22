%% Start
BEGIN_ACADO;
acadoSet('problem_name','compressor_control_2D');
%% Initial setting of the model
DifferentialState x1 x2; Control u;

f = acado.DifferenntialEquation();
f.add(dot(x1) == 1/B*(x2-gamma*((x1+x10)^(1/2)-(x10)^(1/2))));
f.add(dot(x2) == B*(-k3*x2^3-k2*x2^2-k1*x2-x1-u));

%% Initial setting of the optimisation problem
S = 1e2*eye(3); 
S(1,1) = 1e2; 
h = {u x1 x2} ; r = [0 0 0] ;

ocp = acado.OCP(0.0, 20, 50) ;
ocp.minimizeLSQ(S , h , r) ;
ocp.subjectTo(f);
ocp. subjectTo( 'AT START' , x1 == 0.1); 
ocp.subjectTo('AT START ', x2 == 0);
ocp.subjectTo(0 <= u <= 0.2);

%% Initial setup of the acado algorithm
algo=acado.OptimizationAlgorithm(ocp);
algo.set('KKT TOLERANCE', 1e-6);
algo.set('INTEGRATOR TYPE', 'INT RK45');
algo.set('HESSIAN APPROXIMATION', 'GAUSS NEWTON');
algo.set('DISCRETIZATION TYPE', 'COLLOCATION');
END_ACADO;
%% Run
out = compressor_contorl_2D_RUN();