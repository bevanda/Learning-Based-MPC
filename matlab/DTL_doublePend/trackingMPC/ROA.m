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

%% Until N do

% do the first step N=1
F_w=[F_f*A F_f*B;...
    F_x zeros(lenX,m);...
    zeros(lenU,n) F_u];
h_w=[h_f;h_x;h_u];
pre_set = Polyhedron(F_w,h_w);
pre_set=projection(pre_set,1:n);
for i=2:N
    F_p = pre_set.A;
    h_p = pre_set.b;
    F = [F_p*A F_p*B;...
        F_x zeros(lenX,m);...
        zeros(lenU,n) F_u];
    h = [h_p;h_x;h_u];
    pre_set=Polyhedron(F,h);
    pre_set=projection(pre_set,1:n);
end
roa=pre_set;
end