import casadi.*


x = MX.sym('x',2); % Two states
p = MX.sym('p');   % Free parameter

% Expression for ODE right-hand side
z = 1-x(2)^2;
rhs = [z*x(1)-x(2)+2*tanh(p);x(1)];

% ODE declaration with free parameter
ode = struct('x',x,'p',p,'ode',rhs);

% Construct a Function that integrates over 1s
F = integrator('F','cvodes',ode,struct('tf',1));

% Control vector
u = MX.sym('u',5,1);

x = [0;1]; % Initial state
for k=1:5
  % Integrate 1s forward in time:
  % call integrator symbolically
  res = F('x0',x,'p',u(k));
  x = res.xf;
end

% NLP declaration
nlp = struct('x',u,'f',dot(u,u),'g',x);

% Solve using IPOPT
solver = nlpsol('solver','ipopt',nlp);
res = solver('x0',0.2,'lbg',0,'ubg',0);

plot(full(res.x))