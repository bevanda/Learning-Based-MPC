classdef casadi_block < matlab.System & matlab.system.mixin.Propagates
    % untitled Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    properties
        % Public, tunable properties.

    end

    properties (DiscreteState)
    end

    properties (Access = private)
        % Pre-computed constants.
        casadi_solver
        x0
        lbx
        ubx
        lbg
        ubg
    end

    methods (Access = protected)
        function num = getNumInputsImpl(~)
            num = 2;
        end
        function num = getNumOutputsImpl(~)
            num = 1;
        end
        function dt1 = getOutputDataTypeImpl(~)
        	dt1 = 'double';
        end
        function dt1 = getInputDataTypeImpl(~)
        	dt1 = 'double';
        end
        function sz1 = getOutputSizeImpl(~)
        	sz1 = [1,1];
        end
        function sz1 = getInputSizeImpl(~)
        	sz1 = [1,1];
        end
        function cp1 = isInputComplexImpl(~)
        	cp1 = false;
        end
        function cp1 = isOutputComplexImpl(~)
        	cp1 = false;
        end
        function fz1 = isInputFixedSizeImpl(~)
        	fz1 = true;
        end
        function fz1 = isOutputFixedSizeImpl(~)
        	fz1 = true;
        end
        function setupImpl(obj,~,~)
            % Implement tasks that need to be performed only once, 
            % such as pre-computed constants.
            
            import casadi.*

            T = 0.5; % Time horizon
            N = 50; % number of control intervals
            % eqilibilium point
            x_eq = [0.500000000000000;1.68750000000000;1.15470000000000;0];
            u_eq = 1.15470000000000;            
            x_init = [0.150000000000000;1.28750000000000;1.15470000000000;0];

            % Constraints of the compressor model
            mflow_min=0; mflow_max=1;
            prise_min=1.1875; prise_max=2.1875;
            throttle_min=0.1547; throttle_max=2.1547;
            throttle_rate_min=-20; throttle_rate_max=20;
            u_min=0.1547;u_max=2.1547;

            umax = u_max; umin = u_min;
            xmax = [mflow_max; prise_max; throttle_max; throttle_rate_max]; 
            xmin = [mflow_min; prise_min; throttle_min; throttle_rate_min];

            
            % Declare model variables
            x1 = SX.sym('x1');
            x2 = SX.sym('x2');
            x3 = SX.sym('x3');
            x4 = SX.sym('x4');
            x = [x1; x2; x3; x4];
            u = SX.sym('u');

            % Model equations
            xdot = [-x2+1+3*(x1/2)-(x1^3/2);...       % mass flow rate 
                    (x1+1-x3*sqrt(x2));...            % pressure rise rate 
                    x4;...                                % throttle opening rate
                    -1000*x3-2*sqrt(500)*x4+1000*u];    % throttle opening accelerat


            % Objective term
            L = (x1-x_eq(1))^2 + (x2-x_eq(2))^2 + (x3-x_eq(3))^2 + (x4-x_eq(4))^2  + (u-u_eq)^2;
            
            % Continuous time dynamics
            f = casadi.Function('f', {x, u}, {xdot, L});

            % Formulate discrete time dynamics
            % Fixed step Runge-Kutta 4 integrator
            M = 4; % RK4 steps per interval
            DT = T/N/M;
            f = Function('f', {x, u}, {xdot, L});
            X0 = MX.sym('X0', 4);
            U = MX.sym('U');
            X = X0;
            Q = 0;
            for j=1:M
               [k1, k1_q] = f(X, U);
               [k2, k2_q] = f(X + DT/2 * k1, U);
               [k3, k3_q] = f(X + DT/2 * k2, U);
               [k4, k4_q] = f(X + DT * k3, U);
               X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
               Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
            end
            F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'});

            % Start with an empty NLP
            w={};
            w0 = [];
            lbw = [];
            ubw = [];
            J = 0;
            g={};
            lbg = [];
            ubg = [];

            % "Lift" initial conditions
            X0 = MX.sym('X0', 4);
            w = {w{:}, X0};
            lbw = [lbw; x_init];
            ubw = [ubw; x_init];
            w0 = [w0; x_init];

            % Formulate the NLP
            Xk = X0;
            for k=0:N-1
                % New NLP variable for the control
                Uk = MX.sym(['U_' num2str(k)]);
                w = {w{:}, Uk};
                lbw = [lbw; umin];
                ubw = [ubw;  umax];
                w0 = [w0;  0];

                % Integrate till the end of the interval
                Fk = F('x0', Xk, 'p', Uk);
                Xk_end = Fk.xf;
                J=J+Fk.qf;

                % New NLP variable for state at end of interval
                Xk = MX.sym(['X_' num2str(k+1)], 4);
                w = {w{:}, Xk};
                lbw = [lbw; xmin];
                ubw = [ubw; xmax];
                w0 = [w0; zeros(4,1)];

                % Add equality constraint
                g = {g{:}, Xk_end-Xk};
                lbg = [lbg; zeros(4,1)];
                ubg = [ubg; zeros(4,1)];
            end

            % Create an NLP solver
            prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
            options = struct('ipopt',struct('print_level',0),'print_time',false);
            solver = nlpsol('solver', 'ipopt', prob, options);

            obj.casadi_solver = solver;
            obj.x0 = w0;
            obj.lbx = lbw;
            obj.ubx = ubw;
            obj.lbg = lbg;
            obj.ubg = ubg;
        end

        function u = stepImpl(obj,x,t)
            disp(t)
            tic
            w0 = obj.x0;
            lbw = obj.lbx;
            ubw = obj.ubx;
            solver = obj.casadi_solver;
            lbw(1:4) = x;
            ubw(1:4) = x;
            sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
                        'lbg', obj.lbg, 'ubg', obj.ubg);
  
            u = full(sol.x(5));
            toc
        end

        function resetImpl(obj)
            % Initialize discrete-state properties.
        end
    end
end
