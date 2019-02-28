function [x,f] = L2NW(xt,yt,dt,kernel,bandwidth)

% NPREGRESS_NW Performs nonparametric Nadaraya-Watson kernel regression
%   [x,f,pval] = L2NW(xt,yt,kernel,nbootstraps,dt) performs 
% nonparametric regression of 'yt' onto 'xt', smoothed using the kernel 
% specified by 'kernel'.
%
% 'nbootstraps' specifies the number of bootstrap repetitions used to 
% compute standard error of the mean of the estimator.
if nargin<7, nbootstraps = []; end
if nargin<6, nbins = []; end
if nargin<5, bandwidth = []; end
if nargin<4, kernel = []; end
if nargin<3, dt = []; end
compute_sem =0;
if isempty(dt), dt = 1; end
n = length(xt);
if isempty(kernel), kernel = 'Epanechnikov'; end
if isempty(bandwidth), bandwidth = (2.5/100)*range(xt); end % kernel bandwidth is 2.5% of the total range
if isempty(nbins), nbins = round(range(xt)/bandwidth); end % values of x at which to estimate f(x) -- 1 point per bandwidth
if ~isempty(nbootstraps), compute_sem = 1; end

%% define kernel function
if strcmp(kernel,'Uniform'), kernel = @(x,mu,bandwidth) (1/2)*(abs((x - mu)/bandwidth) < 1);
elseif strcmp(kernel,'Epanechnikov'), kernel = @(x,mu,bandwidth) (3/4)*(1-((x - mu)/bandwidth).^2).*(abs((x - mu)/bandwidth) < 1);
elseif strcmp(kernel,'Biweight'), kernel = @(x,mu,bandwidth) (15/16)*((1-((x - mu)/bandwidth).^2).^2).*(abs((x - mu)/bandwidth) < 1);
elseif strcmp(kernel,'Gaussian'), kernel = @(x,mu,bandwidth) exp(-(((x - mu)/bandwidth).^2)/2);
end

%% determine tuning function
if ~compute_sem % just return the means
    binedges = linspace(min(xt),max(xt),nbins+1);
    xval = 0.5*(binedges(1:end-1) + binedges(2:end));
    kval = cell2mat(arrayfun(@(xi) kernel(xi,xt,bandwidth),xval,'UniformOutput',false));
    x.mu = xval;
    f.mu = sum(kval.*repmat(yt,[1 nbins]))./(sum(kval)*dt);
else % obtain both mean and sem by bootstrapping (slow)
    x_mu = zeros(nbootstraps,nbins);
    f_mu = zeros(nbootstraps,nbins);
    for j=1:nbootstraps
        sampindx = sort(randsample(1:n,n,true));  % sample with replacement
        xt_samp = xt(sampindx); yt_samp = yt(sampindx);
        binedges = linspace(min(xt_samp),max(xt_samp),nbins+1);
        xval = 0.5*(binedges(1:end-1) + binedges(2:end));
        kval = cell2mat(arrayfun(@(xi) kernel(xi,xt_samp,bandwidth),xval,'UniformOutput',false));
        x_mu(j,:) = xval;
        f_mu(j,:) = sum(kval.*repmat(yt_samp,[1 nbins]))./(sum(kval)*dt);
    end
    x.mu = mean(x_mu);
    x.sem = std(x_mu);
    f.mu = mean(f_mu); % mean
    f.sem = std(f_mu); % standard error of the mean
end

