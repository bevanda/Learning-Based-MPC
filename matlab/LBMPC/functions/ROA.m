function [Vrep]=ROA(Ak,W,N)
% get H-representations of polyhedra
F_f=W.A; h_f=W.b;
n=size(Ak,1);
%% Doing forward reachability N steps:
%  R = { x | A*x+B*u \in X, u \in U } starting from X_f
F_w=F_f*Ak;
h_w=h_f;
predecessor = Polyhedron(F_w,h_w);
predecessor=projection(predecessor,1:n);
precesessor.minVRep();
for i=2:N
    F_p = predecessor.A;
    h_p = predecessor.b;
    F = F_p*Ak;
    h = h_p;
    predecessor=Polyhedron(F,h);
    predecessor=projection(predecessor,1:n);
end
Vrep=predecessor;
end