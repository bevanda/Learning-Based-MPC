function    [F_x,h_x,...    % nominal state constraints polytope
            F_u,h_u,...     % nominal input constraints polytope
            F_w_N,h_w_N...  % terminal extended state constraints polytope
            ]...
            =getCONS(...
                xmax,xmin,umax,umin,...
                x_wp,u_wp,m,n,...
                A,B,K,LAMBDA,PSI,LAMBDA_0,PSI_0...
                )
%==========================================================================
% Defining polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x (so called H-representation)
%==========================================================================
disp('Computing constraint polytopes...');
% Shift the system constraints w.r.t. to the linearisation point
F_u = [eye(m); -eye(m)]; h_u = [umax-u_wp; -umin+u_wp];
F_x = [eye(n); -eye(n)]; h_x = [xmax-x_wp; -xmin+x_wp];

% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);

%==========================================================================
% Compute maximal positively invariant set (MPIS)
%==========================================================================
lambda=0.99; % can be chosen arbitrarily close to 1, the obtained
% invariant set can be used as a reliable polyhedral approximation of MPIS 
disp('Computing and simplifying terminal set...');
% tic;
% extended state constraints
% L=(PSI - K*LAMBDA);
% L0=(PSI_0-K*LAMBDA_0);
% F_w = [F_x zeros(length_Fx, m);
%     zeros(length_Fx, n) F_x*LAMBDA; ...
%     F_u*K, F_u*L; ...
%     zeros(length_Fu, n) F_u*PSI];
% h_w = [...
%     h_x; ...
%     lambda*(h_x - F_x*LAMBDA_0); ...
%     h_u-F_u*L0; ...
%     lambda*(h_u- F_u*PSI_0)];
% 
% Xw = Polyhedron(F_w, h_w);
% Aw=[A+B*K B*L; zeros(m,n) eye(m)];
% 
% X_w_N=compute_MPIS(Xw,Aw);
% 
% disp('Terminal set Polyhedron:');
% X_w_N %#ok
% X_w_N.minHRep(); % simplifying the polyhedron/constraints
% F_w_N = X_w_N.A; h_w_N = X_w_N.b;
% toc
set=load('../saved_data+plots/data/term_set.mat');
F_w_N = set.F_w_N ; h_w_N = set.h_w_N;
% term_poly=Polyhedron(F_w_N,h_w_N);
% term_poly.projection(1:2);
% term_poly.plot('wire',1,'linewidth',2.5,'linestyle',':','color', 'blue');
end