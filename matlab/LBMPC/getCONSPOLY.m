function    [F_x,h_x,...    % nominal state constraints polytope
            F_u,h_u,...     % nominal input constraints polytope
            F_w_N,h_w_N...  % terminal extended state constraints polytope
            ,F_x_d,h_x_d]... % uncertainty polytope
            =getCONSPOLY(xmax,xmin,umax,umin,state_uncert,x_wp,u_wp,m,n,A,B,Q,R,LAMBDA,PSI,LAMBDA_0,PSI_0)
%==========================================================================
% Define polytopic constraints on input F_u*x <= h_u and
% state F_x*x <= h_x and define model uncertainty as a F_g*x <= h_g
% (so called H-representation)
%==========================================================================
% Shift the system constraints w.r.t. to the linearisation point
F_u = [eye(m); -eye(m)]; h_u = [umax-u_wp; -umin+u_wp];
F_x = [eye(n); -eye(n)]; h_x = [xmax-x_wp; -xmin+x_wp];
F_g = [eye(n); -eye(n)]; h_g = [state_uncert; state_uncert]; % uncertainty polytope
% count the length of the constraints on input, states, and uncertainty:
length_Fu = length(h_u);
length_Fx = length(h_x);
length_Fg = length(h_g);

% uncertainty polytope
temp = Polyhedron(F_x, h_x) - Polyhedron(F_g, h_g);
temp.minHRep();
F_x_d= temp.A; h_x_d = temp.b;
%==========================================================================
% Compute maximal positively invariant set (MPIS)
%==========================================================================

% Terminal feedback policy for terminal set computations
% r_i as the inverse of the square of the maximum permissible value 
% for the corresponding u_i
maxadm_controlweight = 10; 
K_t = -dlqr(A, B, Q, maxadm_controlweight*R);
lambda=0.99; % can be chosen arbitrarily close to 1, the obtained
% invariant set can be used as a reliable polyhedral approximation of MPIS 
disp('Computing and simplifying terminal set...');
% extended state constraints
L=(PSI - K_t*LAMBDA);
L0=(PSI_0 - K_t*LAMBDA_0);
F_w = [ F_x zeros(length_Fx, m);
        zeros(length_Fx, n) F_x*LAMBDA; ...
        F_u*K_t, F_u*L; ...
        zeros(length_Fu, n) F_u*PSI; ...
        F_x_d*(A+B*K_t) F_x_d*B*L];
h_w = [ h_x; ...
        lambda*(h_x - F_x*LAMBDA_0); ...
        h_u - F_u*L0; ...
        lambda*(h_u - F_u*PSI_0); ...
        h_x_d - F_x_d*B*(PSI_0-K_t*LAMBDA_0)];

% disturbance constraints of the extended state 
F_g_w = [F_g zeros(length_Fg,m); ...
        zeros(m, n) eye(m); ...
        zeros(m, n) -eye(m)];
h_g_w = [h_g; ...
        zeros(2*m,1)];
    
% calculating the robust positively invariant set    
[F_w_N0, h_w_N0] = pdiff(F_w, h_w, F_g_w, h_g_w);

term_poly = Polyhedron(F_w_N0, h_w_N0);
term_poly.minHRep(); % simplifying the polyhedron/constraints
F_w_N = term_poly.A; h_w_N = term_poly.b;
disp('Terminal set Polyhedron:');
term_poly %#ok

end