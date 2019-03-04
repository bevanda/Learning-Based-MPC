function [g] = L2NWoracle(ksi,X,Y,kernel,bandwidth,lambda)

% [g] = L2NW(xi,kernel,bandwidth,lambda)  performs ...
% nonparametric L2 normalised Nadaraya-Watson kernel regression
% xi = [x' u']' 
% 'lambda' term ensures regularisation and differentiabilty
% 'bandwidith' term sets up the smoothnes of the estimator'
if nargin < 4, kernel = []; end
if nargin < 5, bandwidth = []; end
if nargin < 6, lambda = []; end
if isempty(kernel), kernel = 'Epanechnikov'; end
if isempty(bandwidth), bandwidth = 0.5; end 
if isempty(lambda),lambda = 0.001; end
%% Defeining CITE 3D kernel

if strcmp(kernel,'Epanechnikov'), kernel = @(x,mu,bandwidth) (3/4)*(1-((x - mu)/bandwidth).^2).*(abs((x - mu)/bandwidth) );
end
%% Tuning function setup
kval = arrayfun(@(xi) kernel(xi,ksi,bandwidth),X,'UniformOutput', false);
cell2mat(kval(1))
% 'repmat(yt,[1 nbins])' stacking yt for every data point
g = sum(kval.*Y)./(lambda+sum(kval));
end