LAMBDA = [0.342368418956904;-0.385164471326517;0.395332936252230;-2.22651758141623e-15];
PSI = [0.395332936252230];
x_s = [0];
f = @(x,theta) (x-LAMBDA*theta).'*eye(4)*(x-LAMBDA*theta)+(x_s-LAMBDA*theta).'*(1000*eye(4))*(x_s-LAMBDA*theta);
fsurf(f,'ShowContours','on');