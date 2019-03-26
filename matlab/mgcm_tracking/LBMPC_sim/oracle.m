function g = oracle(ksi, X, Y)
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
for i=1:n
    % Gaussian Kernel
    kval(i) = exp(-((norm(X(:,i) - ksi).^2)/bandwidth));
end
skval = sum(kval);
weight = kval/(lambda + skval);
y=zeros(s,1);
for i=1:n
    out = weight(i)*Y(:,i);
    y = y + out;
end
g = y;
