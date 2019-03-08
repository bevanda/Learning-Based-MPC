function [alpha, lambda] = mpcalpha_cn(c,N,varargin)
% mpcalpha_cn([c0, c1, ...], N, m, endweight, output)
% Computes the suboptimality index alpha for an MPC problem
% satisfying the controllability condition with 
% beta(r,n) = cn*r for optimization horizon N
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
% (C) Lars Gruene 2007

if (nargin < 2)
    fprintf('mpcalpha_cn needs at least two arguments\n');
    fprintf('see "help mpcalpha_cn" for more information\n');
    return;
end;
    
if (nargin>=3)
    m = varargin{1};
else
    m = 1;
end;

if (nargin>=4) 
    endweight = varargin{2};
else
    endweight = 1;
end;

if (nargin>=5)
    output = varargin{3};
else
    output = 0;
end;

% extend c sequence, if shorter than N
c = [c, zeros(1,N-length(c))];

% compute B_k
g = cumsum(c) + (endweight-1.0)*c;

% inequality constraint matrix
A = zeros(2*N-m-1,N+1);
for i=1:N-1
  A(i,i) = 1-g(N-i+1);
  for j=i+1:N-1
    A(i,j) = 1;
  end
  A(i,N) = endweight;
end;
for i=m+1:N
  for j=m+1:i-1 
     A(N+i-m-1,j) = -1; 
  end;
  A(N+i-m-1,i) = A(N+i-m-1,i)-g(N-i+m+1);
  A(N+i-m-1,N+1) = 1;  
end;

if output>=2
    A
end;

% right hand side of inequalities
b = zeros(2*N-m-1,1);

% lower bound 0, no upper bound
lb = zeros(N+1,1);
ub = [];

% objective vector
f = [ones(N-1,1);endweight;-1];

% equality constraints
Aeq = [ ones(1,m), zeros(1,N-m+1)];
beq = [1];

if (output>=2)
%   options = optimset('LargeScale', 'on', 'Simplex', 'off'); % Simplex
    options = optimset('display', 'off','LargeScale', 'on');
%   dpoes not exis anymore
else
  options = optimset('display', 'off','LargeScale', 'on');
end;
  
% solve linear program, l = [lambda_0,...,lambda_{N-1},nu]
l = linprog(f,A,b,Aeq,beq,lb,ub,[],options);

% compute alpha and value functions
alpha = f'*l;

VNp  = l(N+1);
VN   = sum(l(1:N-1))+endweight*l(N);
lambda = l(1:N);  

% plot worst case sequences
if (output>=1)
  ep = N-m;
  
  lp = zeros(N,1);
  for j=1:N
      if j<ep
          lp(j) = l(j+m);
      else
          lp(j) = c(j-ep+1)*l(ep+m);
      end;
  end;
  
%  l(1:N)
  
  pl = [l(1:N);lp];

  plot([0:N-1],l(1:N),'.',[m:N+m-1],lp,'ro');
  axis([-0.1,N+m-0.9,-0.1,max(max(pl),1)*1.1]);
end;

if (output>=2)
  fprintf('V(%d,x(0)) = %f\n',N,VN);
  fprintf('V(%d,x(%d)) = %f\n',N,m,VNp);

  fprintf('alpha     = %f\n', alpha)
  fprintf('OptimBound= %f\n', 1/alpha);

  if (alpha<=0) 
    fprintf('UNSTABLE\n');
  else
    fprintf('STABLE\n');
  end;
end;

