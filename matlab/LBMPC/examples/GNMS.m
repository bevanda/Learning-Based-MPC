classdef GNMS < handle
  properties
    % Symbolic representation of the OCP
    ocp
    % OCP data
    data
    % Solver options
    opts
    % Verbose output
    verbose
    % Dimensions
    N
    nx
    nu
    % Solution time grid
    t
    % Control interval length
    dt
    % CasADi function objects
    fun
    % Nonlinear program
    nlp
    % Current iteration
    n_iter
    % Current solution
    sol
  end
  methods
    function self = GNMS(ocp, data, opts)
      % Constructor
      self.ocp = ocp;
      self.data = data;
      self.opts = opts;

      % Verbosity?
      if (~isfield(opts, 'verbose'))
	self.verbose = opts.verbose;
      else
	self.verbose = true;
      end

      % Problem dimensions
      self.N = self.opts.N;
      self.nx = numel(self.ocp.x);
      self.nu = numel(self.ocp.u);

      % Time grid
      self.t = linspace(0, data.T, self.N+1);

      % Interval length
      self.dt = data.T / self.N;

      % Get discrete-time dynamics
      self.rk4();

      % Transcribe to NLP
      self.transcribe();

      % Initialize x
      self.sol.w = self.nlp.w0;

      % Extract x,u from w
      [x_traj,u_traj] = self.fun.traj(self.sol.w);
      self.sol.x = full(x_traj);
      self.sol.u = full(u_traj);

      % Initialize Gauss-Newton method
      self.init_gauss_newton();
    end

    function init_gauss_newton(self)
      % Linearize the problem w.r.t. w
      self.nlp.J = jacobian(self.nlp.g, self.nlp.w);
      self.nlp.JM = jacobian(self.nlp.M, self.nlp.w);

      % Gauss-Newton Hessian and gradient
      self.nlp.H = self.nlp.JM' * self.nlp.JM;
      self.nlp.c = self.nlp.JM' * self.nlp.M;

      % Display sparsities
      disp('J sparsity pattern:')
      self.nlp.J.sparsity().spy();
      disp('H sparsity pattern:')
      self.nlp.H.sparsity().spy();

      % Functions for calculating g, J, H
      self.fun.g = casadi.Function('g', {self.nlp.w}, {self.nlp.g},...
                                   {'w'}, {'g'});
      self.fun.J = casadi.Function('J', {self.nlp.w}, {self.nlp.J},...
                                   {'w'}, {'J'});
      self.fun.H = casadi.Function('H', {self.nlp.w}, {self.nlp.H, self.nlp.c},...
                                   {'w'}, {'H', 'c'});

      % Out new step is determined by solving the quadratic program for
      % dw = w_new-w
      %        minimize    1/2 dw'*H*dw + c'*dw
      %        subject to  g + A*dw = 0
      %                    lbw-w <= dw <= ubw-w
      qp = struct('h', self.nlp.H.sparsity(), 'a', self.nlp.J.sparsity());
      qp_options = struct();
      if ~self.verbose
        qp_options.printLevel = 'none';
      end
      self.fun.qp_solver = casadi.conic('qp_solver', 'hpmpc', qp, qp_options);

      % Iteration counter
      self.n_iter = 0;

      % Residual
      self.sol.norm_dw = inf;
    end

    function rk4(self)
        % Continuous-time dynamics
        x = self.ocp.x;
        u = self.ocp.u;
        ode = self.ocp.ode;
        f = casadi.Function('f', {x, u}, {ode}, {'x','p'}, {'ode'});

        % Implement RK4 integrator that takes a single step
        k1 = f(x, u);
        k2 = f(x+0.5*self.dt*k1, u);
        k3 = f(x+0.5*self.dt*k2, u);
        k4 = f(x+self.dt*k3, u);
        xk = x+self.dt/6.0*(k1+2*k2+2*k3+k4);

        % Return as a function
        self.fun.F = casadi.Function('RK4', {x,u}, {xk}, {'x0','p'}, {'xf'});

	% Least squares objective function
	lsq = self.ocp.lsq;
	self.fun.Lsq = casadi.Function('Lsq', {x, u}, {lsq}, {'x','p'}, {'lsq'});
    end

    function transcribe(self)
      % Start with an empty NLP
      w = {}; % Variables
      g = {}; % Equality constraints
      M = {}; % Least squares objective function
      lbw = {}; % Lower bound on w
      ubw = {}; % Upper bound on w
      w0 = {}; % Initial guess for w

      % Expressions corresponding to the trajectories we want to plot
      x_plot = {};
      u_plot = {};

      % Initial conditions
      xk = casadi.MX.sym('x0', self.nx);
      w{end+1} = xk;
      lbw{end+1} = self.data.x0;
      ubw{end+1} = self.data.x0;
      w0{end+1} = self.data.x_guess;
      x_plot{end+1} = xk;

      % Loop over all times
      for k=0:self.N-1
          % Define local control
          uk = casadi.MX.sym(['u' num2str(k)], self.nu);
          w{end+1} = uk;
          lbw{end+1} = self.data.u_min;
          ubw{end+1} = self.data.u_max;
          w0{end+1} = self.data.u_guess;
          u_plot{end+1} = uk;

          % Simulate the system forward in time
          Fk = self.fun.F('x0', xk, 'p', uk);
    	  x_next = Fk.xf;

          % Add least squares term to the objective
          M{end+1} = self.fun.Lsq(xk, uk);

          % Define state at the end of the interval
          xk = casadi.MX.sym(['x' num2str(k+1)], self.nx);
          w{end+1} = xk;
          if k==self.N-1
              lbw{end+1} = self.data.xN;
              ubw{end+1} = self.data.xN;
          else
              lbw{end+1} = self.data.x_min;
              ubw{end+1} = self.data.x_max;
          end
          w0{end+1} = self.data.x_guess;
          x_plot{end+1} = xk;

          % Impose continuity
          g{end+1} = xk - x_next;
      end

      % Concatenate variables and constraints
      self.nlp.w = vertcat(w{:});
      self.nlp.g = vertcat(g{:});
      self.nlp.M = vertcat(M{:});
      self.nlp.lbw = vertcat(lbw{:});
      self.nlp.ubw = vertcat(ubw{:});
      self.nlp.w0 = vertcat(w0{:});

      % Create a function that maps the NLP decision variable to the x and u trajectories
      self.fun.traj = casadi.Function('traj', {self.nlp.w}, {horzcat(x_plot{:}), horzcat(u_plot{:})},...
	                             {'w'}, {'x', 'u'});
    end

    function sqpstep(self)
      % Update iteration counter
      self.n_iter = self.n_iter + 1;

      % Calculate the QP matrices
      self.sol.g = self.fun.g(self.sol.w);
      self.sol.J = self.fun.J(self.sol.w);
      [self.sol.H, self.sol.c] = self.fun.H(self.sol.w);

      % Solve the QP to get the step in in w
      qp_solution = self.fun.qp_solver('a', self.sol.J, 'h', self.sol.H, 'g', self.sol.c,...
				       'lbx', self.nlp.lbw-self.sol.w,...
				       'ubx', self.nlp.ubw-self.sol.w,...
				       'lba', -self.sol.g, 'uba', -self.sol.g, 'x0', 0);
      dw = full(qp_solution.x);

      % Check convergence criteria
      self.sol.norm_dw = norm(dw);

      % Take (full) step
      self.sol.w = self.sol.w + dw;

      % Extract x,u from w
      [x_traj,u_traj] = self.fun.traj(self.sol.w);
      self.sol.x = full(x_traj);
      self.sol.u = full(u_traj);

      % Print progress
      if self.verbose || mod(self.n_iter,5)==1
        disp(repmat('-', 1, 70))
        fprintf('%15s %15s\n', 'SQP iteration', 'norm(dw)');
      end
      fprintf('%15d %15g\n', self.n_iter, self.sol.norm_dw);
    end
  end
end