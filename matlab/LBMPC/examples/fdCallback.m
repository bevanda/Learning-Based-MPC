
classdef fdCallback < casadi.Callback
  properties
    d
  end
  methods
    function self = fdCallback(name, d, options)
      self@casadi.Callback();
      self.d = d;
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
      f = sin(self.d * x);
      arg = {f};
    end
  end
end
