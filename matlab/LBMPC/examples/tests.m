%% test function
% Use the function
import casadi.*

f = fdCallback('f', 0.5, struct);
res = f(2);
disp(res)
% options = struct;
% % options.common_options.helper_options = struct('enable_fd',true,'fd_method','forward');
% options.fd_options = struct('enable_fd',true,'fd_method','forward');
% fd_options = struct('enable_fd',true,'fd_method','forward');

x = MX.sym('x');
u = MX.sym('u');
disp(f(x));

% Derivates OPTION 1: finite-differences
eps = 1e-5;
disp((f(2+eps)-f(2))/eps);

f = fdCallback('f', 0.5, struct('enable_fd',true));
J = Function('J',{x},{jacobian(f(x),x)});
disp(J(2));

eps = 1e-5;
disp((f(2+eps)-f(2))/eps);
% o = oracleCallback('o', data, struct('enable_fd',true));

x = SX.sym('x');
u = SX.sym('u');

k{1}=x;
k{2}=u;
k{3}=x^2;
k{4}=u^2;
k{5}=x^3;
sum = 0;
for i=1:5
    sum = sum +k{i};
end
