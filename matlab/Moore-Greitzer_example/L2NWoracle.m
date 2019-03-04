function [g] = L2NWoracle(ksi,X,Y,kernel,bandwidth,lambda)

% [g] = L2NW(xi,kernel,bandwidth,lambda)  performs ...
% nonparametric L2 normalised Nadaraya-Watson kernel regression
% xi = [x' u']' 
% 'lambda' term ensures regularisation and differentiabilty
% 'bandwidith' term sets up the smoothnes of the estimator'
if nargin < 4, kernel = []; end
if nargin < 5, bandwidth = []; end
if nargin < 6, lambda = []; end
if isempty(kernel), kernel = 'Gaussian'; end
if isempty(bandwidth), bandwidth = 0.5; end 
if isempty(lambda),lambda = 0.001; end
d = numel(ksi);
n = size(X,2);
%% Defeining CITE 3D kernel

if strcmp(kernel,'Epanechnikov'), kernel = @(x,mu,bandwidth) (3/4)^d*(1-((x - mu)/bandwidth).^2).*(abs((x - mu)/bandwidth));
elseif strcmp(kernel,'Gaussian'), kernel = @(x,mu,bandwidth) exp(-((norm(x - mu).^2)/bandwidth));
end
%% Tuning function setup

kval = zeros(1,n);
for i=1:n
    kval(i) = exp(-((norm(X(:,i) - ksi).^2)/bandwidth));
end
skval = sum(kval);
weight = kval/(lambda + skval);
g=0;
for i=1:n
    out = weight(i)*Y(:,i);
    g = g + out;
end
end