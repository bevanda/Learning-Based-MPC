
%%
rng(0,'twister'); % For reproducibility
n = 1000;
x = [linspace(-10,10,n)',linspace(-10,10,n)'];
y = 1 + x*5e-2 + sin(x)./x + 0.2*randn(n,2);

gprMdl = fitrgp(x,y,'Basis','linear',...
      'FitMethod','exact','PredictMethod','exact');
  
ypred = resubPredict(gprMdl);

plot(x,y,'b.');
hold on;
plot(x',ypred,'r','LineWidth',1.5);
xlabel('x');
ylabel('y');
legend('Data','GPR predictions');
hold off;

%%
function [X,Y,h] = kernelRegression(x, y, edges, h)
% 
% [X,Y,h] = kernelRegression(x, y, edges, h)
%
% KSR   Kernel smoothing regression
% r=ksr(x,y) returns the Gaussian kernel regression in structure r such that
%   r.f(r.x) = y(x) + e
%
% Algorithm
% The kernel regression is a non-parametric approach to estimate the
% conditional expectation of a random variable:
%
% E(Y|X) = f(X)
%
% where f is a non-parametric function. Based on the kernel density
% estimation, this code implements the Nadaraya-Watson kernel regression
% using the Gaussian kernel as follows:
%
% f(x) = sum(kerf((x-X)/h).*Y)/sum(kerf((x-X)/h))
%
% See also gkde, ksdensity
%
% Adapted from the ksr function y Yi Cao at Cranfield University.

x=x(:);
y=y(:);

% clean missing or invalid data points (nan are removed).

inv=(x~=x)|(y~=y);
x(inv)=[];
y(inv)=[];

% Default parameters
if ~exist('edges','var')
    disp('Using default edges...')
end

N = length(edges);
if ~exist('h','var')
    % optimal bandwidth suggested by Bowman and Azzalini (1997) p.31
    hx=median(abs(x-median(x)))/0.6745*(4/3/length(x))^0.2;
    hy=median(abs(y-median(y)))/0.6745*(4/3/length(x))^0.2;
    h=sqrt(hy*hx);
    if h<sqrt(eps)*N
        error('There is no enough variation in the data. Regression is meaningless.')
    end
    fprintf(2,'Using optimal bandwidth (Bowman and Azzalini, 1997). h = %f.\n',h);
elseif ~isscalar(h)
    error('h must be a scalar.')
end

% Gaussian kernel function
KERNEL = 'gaussian';
% TODO: Implement other kernels i.e. uniform, triangular, epanechnickov.
if strcmp(KERNEL,'gaussian')
    kerf=@(z)exp(-z.*z/2)/sqrt(2*pi);
end

X = edges;
Y = nan(1,N);
for k=1:N
    z=kerf((X(k)-x)/h);
    Y(k)=sum(z.*y)/sum(z);
end
end