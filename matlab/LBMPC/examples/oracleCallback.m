
classdef oracleCallback < casadi.Callback
  properties
  end
  methods
    function self = oracleCallback(name, options)
      self@casadi.Callback();
      construct(self, name, options);
    end

    % Number of inputs and outputs
    function v=get_n_in(self)
      v=1;
    end
    function v=get_n_out(self)
      v=1;
    end

    % Initialize the object
    function init(self)
      disp('initializing object')
    end

    % Evaluate numerically
    function arg = eval(self, arg)
      x = arg{1};
      u = arg{2};
      d = arg{3};
      f = A*x+B*u + oracleL2NW(x,u,d);
      arg = {f};
    end
  end
end
