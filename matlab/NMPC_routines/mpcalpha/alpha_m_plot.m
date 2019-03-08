function alpha_m_plot(C,sigma,N,m,varargin)
% alpha_m_plot(C, sigma, N, m, endweight, subroutineoutput)
% Plots the suboptimality index alpha for an MPC problem
% satisfying the controllability condition with
% beta(r,n) = C*sigma^n*r for optimization horizon N
% for varying control horizons
%
% Arguments:
% C:                   overshoot
% sigma:               decay rate
% N:                   optimization horizon
% m:                   control horizon, default 1
%
% Additional arguments:
% endweight:           weights of last summand in MPC
%                      functional, default 1
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
if (nargin>=5) 
    endweight = varargin{1};
else
    endweight = 1;
end;
if (nargin>=6) 
    subroutineoutput = varargin{2};
else
    subroutineoutput = 0;
end;

% compute alpha values
for i=1:length(m)
    alpha(i) = mpcalpha_Csigma(C,sigma,N,m(i),endweight,subroutineoutput);
end;

figure();
	plot(m,alpha,'*k');
	hold on;
	plot([0,m(length(m))+1],[0 0],'--k');
	axis tight;
	xlabel('Control horizon m');
	ylabel('Suboptimality index \alpha');
	hold off;