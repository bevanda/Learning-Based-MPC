function roa=ROA(sysdata,X_f,Xc,Uc,N)
%% Region of attraciton calculation
% X_f,X_ext: Terminal set, state and input constraints polyhedron

A=sysdata.A;
B=sysdata.B;
n=sysdata.n;
m=sysdata.m;
lenX=length(Xc.b);
lenU=length(Uc.b);
% get H-representations of polyhedra
F_f=X_f.A; h_f=X_f.b;
F_x=Xc.A; h_x=Xc.b;
F_u=Uc.A; h_u=Uc.b;

%% Doing backward reachability N steps:
%  R = { x | A*x+B*u \in X, u \in U } starting from X_f
%
F_w=[F_f*A F_f*B;...
    F_x zeros(lenX,m);...
    zeros(lenU,n) F_u];
h_w=[h_f;h_x;h_u];
predecessor = Polyhedron(F_w,h_w);
predecessor=projection(predecessor,1:n);
for i=2:N
    F_p = predecessor.A;
    h_p = predecessor.b;
    F = [F_p*A F_p*B;...
        F_x zeros(lenX,m);...
        zeros(lenU,n) F_u];
    h = [h_p;h_x;h_u];
    predecessor=Polyhedron(F,h);
    predecessor=projection(predecessor,1:n);
end
roa=predecessor;
end