function [alpha, lambda] = mpcalpha_Csigma(C,sigma,N,varargin)
% mpcalpha_Csigma(C, sigma, N, m, endweight, output)
% Computes the suboptimality index alpha for an MPC problem
% satisfying the controllability condition with 
% beta(r,n) = C*sigma^n*r for optimization horizon N
%
% optional arguments:
% m:         control horizon, default 1
% endweight: weight of last summand in MPC functional, default 1
% output:    =0: no output (default)
%            =1: graphical output
%            =2: graphical and text output
% 
% see the PDF file on 
% www.math.uni-bayreuth.de/~lgruene/publ/mpcbound.html
% for a description of the method, the definition of the 
% controllability condition and the meaning of
% the computed index alpha
%
% The procedure needs "mpcalpha_cn.m" available from
% www.math.uni-bayreuth.de/~lgruene/publ/mpcbound.html
%
% (C) Lars Gruene 2007

if (nargin < 3)
    fprintf('mpcalpha_Csigma needs at least three arguments\n');
    fprintf('see "help mpcalpha_Csigma" for more information\n');
    return;
end;
    

c = C*cumprod([1, sigma*ones(1,N-1)]);

[alpha, lambda] = mpcalpha_cn(c,N,varargin{:});
