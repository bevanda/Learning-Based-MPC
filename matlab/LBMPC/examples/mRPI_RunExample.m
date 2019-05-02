clearvars;

A=[1.01126321746508 -0.0100340214950357 6.46038913508018e-05 1.93716902346107e-07;0.0100340214950357 0.995515380253533 -0.0127681799951143 -5.57226765949307e-05;0 0 0.957038195891878 0.00792982548734093;0 0 -7.92982548734093 0.602405619103784];
K=[-3.07418713694111 2.09578024408840 0.119436236659413 -0.00894688869207666];
B=[4.95338239742909e-07;-0.000193159646826654;0.0429618041081220;7.92982548734093];
Ak=A+B*K;
n=size(A,1);

state_uncert = [0.02;5e-04;0;0];
F_d = [eye(n); -eye(n)]; h_d = [state_uncert; state_uncert];

s_0=10;
eps=5e-5;

polyH = Polyhedron(F_d,h_d); 
polyH.minHRep();
F_w=polyH.A; h_w=polyH.b;

polyV=Polyhedron(polyH.V);
polyR=reach_set(Ak,polyV,s_0);
figure;
polyR=polyR.projection(1:2);
polyR.plot();

mRPIS=calc_mRPIS(Ak,F_w,h_w,eps);
figure;
mRPIS.projecion(1:2);
mRPIS.plot();