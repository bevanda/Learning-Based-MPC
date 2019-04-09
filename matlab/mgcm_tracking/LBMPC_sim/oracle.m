function g = oracle(x,u,data)
X=data.X; Y=data.Y; ksi=[x(1:2);u];
% g = oracle(ksi, X, Y) performs a nonparametric L2 normalised Nadaraya-Watson kernel regression
if nargin <2
    X=zeros(3,1);
    Y=zeros(4,1);
end 

bandwidth = 0.5;  
lambda = 0.001; 
% d = numel(ksi);
s = size(Y,1);
n = size(X,2);

kval = zeros(1,n);
U = @(X,ksi,bandwidth) (norm(X - ksi).^2)/bandwidth^2;

% Gaussian Kernel
gauss = @(X,ksi,bandwidth) exp(-(U(X,ksi,bandwidth)));
% % Spherical Epanechnikov kernel
% epan_sph= @(X,ksi,bandwidth) (1-(U(X,ksi,bandwidth))*(U(X,ksi,bandwidth)'*U(X,ksi,bandwidth) < 1));

for i=1:n
    kval(i) = gauss(X(:,i),ksi,bandwidth);
end
skval = sum(kval);
weight = kval/(lambda + skval);
y=zeros(s,1);
for i=1:n
    out = Y(:,i)*weight(i);
    y = y + out;
end
g = y;
