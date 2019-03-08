function alpha_omega_plot(C,sigma,N,m,endweight,varargin)
% alpha_omega_plot(C, sigma, N, m, endweight, subroutineoutput)
% Plots the suboptimality index alpha for an MPC problem
% satisfying the controllability condition with
% beta(r,n) = C*sigma^n*r for optimization horizon N
% for varying endweights
%
% Arguments:
% C:                   overshoot
% sigma:               decay rate
% N:                   optimization horizon
% m:                   control horizon, default 1
% endweight:           array of weights of last summand
%                      in MPC functional, default 1
%
% Additional arguments:
% subroutineoutput:    = 0: no output (default)
%                      = 1: graphical output
%                      = 2: graphical and text output
%
% see the PDF files on
% www.math.uni-bayreuth.de/~lgruene/publ/mpcbound.html
% www.math.uni-bayreuth.de/~lgruene/publ/mpcvary.html
% for a description of the method, the definition of the
% controllability condition and the meaning of
% the computed index alpha
%
% The procedure needs "mpcalpha_cn.m" and "mpcalpha_Csigma.m"
% available from
% www.math.uni-bayreuth.de/~lgruene/publ/mpcbound.html
%
% (C) Lars Gruene, Juergen Pannek 2010
if (nargin>=6) 
    subroutineoutput = varargin{1};
else
    subroutineoutput = 0;
end;

% compute alpha values
for i=1:length(endweight)
    alpha(i) = mpcalpha_Csigma(C,sigma,N,m,endweight(i),subroutineoutput);
end;

figure();
	plot(endweight,alpha,'*k');
	hold on;
	plot([0,endweight(length(endweight))],[0 0],'--k');
	axis tight;
	xlabel('Endweight \omega');
	ylabel('Suboptimality index \alpha');
	hold off;