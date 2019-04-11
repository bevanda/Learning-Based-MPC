function Theta = SA(Theta_min,Theta_max,Dr,A,B,K_gain,Phaii,Saii,xs,Oracle,R,Q,T,P)
N=size(Dr,2); % Horizon Length
dt=0.001;
iter=0;
D=size(Theta_min,1)*N; % Dimesnion of the problem
% SA Controlling parameters
tempi=1; % Initial Temperature
temps=0.001; % Final Temperature
cf=0.9996; % Cooling Factor
temp=tempi; % Current Temperature
SigmAaa=0.1*(Theta_max(1)-Theta_min(1)); % A parameter in the Gaussian proposal
% Intial solution
Xsol=Theta_min(1)+(rand(1,D)*(Theta_max(1)-Theta_min(1)));
% Point estimate of the target distribution
Theta(1,1:N)=Xsol(1:N); Theta(2,1:N)=Xsol(N+1:D);
U_stable=(Phaii-(-K_gain*Saii))*Theta; U(:,1)=(-K_gain*xs)+U_stable(:,1);
Xsystem(:,1)=xs;
for ij=1:1:N-1
Xsystem(:,ij+1)=( dt*( (A*Xsystem(:,ij)) + (B*U(:,ij)) + (Dr(:,ij)) + Oracle(:,ij) ) ) + Xsystem(:,ij);
U(:,ij+1)=(-K_gain*Xsystem(:,ij+1))+U_stable(:,ij+1);
end
SumX=0; SumU=0;
for ij=2:1:N-1
SumX=SumX+((Xsystem(:,ij)-(Saii*Theta(:,ij)))'*Q*(Xsystem(:,ij)-(Saii*Theta(:,ij))));
SumU=SumU+((U(:,ij)-(Phaii*Theta(:,ij)))'*R*(U(:,ij)-(Phaii*Theta(:,ij))));
end
h_X= ((Xsystem(:,N)-(Saii*Theta(:,N)))'*P*(Xsystem(:,N)-(Saii*Theta(:,N))))...
+((-(Saii*Theta(:,N)))'*T*(-(Saii*Theta(:,N))))+SumX+SumU;
Pai_X=exp(-h_X/temp);
Fmin=h_X;
while temp > temps
iter = iter+1; temp=cf*temp;
%-- Proposal distribution "P(X,Y)" (Gaussian pdf)
Ysol = normrnd(Xsol,SigmAaa);
ind=find(Ysol>Theta_max(1)); Ysol(ind)=Theta_max(1); ind=find(Ysol<Theta_min(1)); Ysol(ind)=Theta_min(1);
% Point estimate of the target distribution
Theta(1,1:N)=Ysol(1:N); Theta(2,1:N)=Ysol(N+1:D);
U_stable=(Phaii-(-K_gain*Saii))*Theta; U(:,1)=(-K_gain*xs)+U_stable(:,1);
Xsystem(:,1)=xs;
for ij=1:1:N-1
Xsystem(:,ij+1)=( dt*( (A*Xsystem(:,ij)) + (B*U(:,ij)) + (Dr(:,ij)) + Oracle(:,ij) ) ) + Xsystem(:,ij);
U(:,ij+1)=(-K_gain*Xsystem(:,ij+1))+U_stable(:,ij+1);
end
SumX=0; SumU=0;
for ij=2:1:N-1
SumX=SumX+((Xsystem(:,ij)-(Saii*Theta(:,ij)))'*Q*(Xsystem(:,ij)-(Saii*Theta(:,ij))));
SumU=SumU+((U(:,ij)-(Phaii*Theta(:,ij)))'*R*(U(:,ij)-(Phaii*Theta(:,ij))));
end
h_Y= ((Xsystem(:,N)-(Saii*Theta(:,N)))'*P*(Xsystem(:,N)-(Saii*Theta(:,N))))...
+((-(Saii*Theta(:,N)))'*T*(-(Saii*Theta(:,N))))+SumX+SumU;
Pai_Y=exp(-h_Y/temp);
%-- Acceptance-Rejection determination
% 1: Draw "u" from Uniform Distribution [0,1]
u=rand;
% 2: Verify to accept or not
alphA=min(1,(Pai_Y/Pai_X));
if u<=alphA % Accept the new state
Xsol=Ysol;
h_X=h_Y;
Pai_X=Pai_Y;
else % Reject
end
% Archive the best solution found so-far
Fmin=min(Fmin, h_X);
Gbest(iter)=Fmin;
% Variation of objective value through the movement of the chain
ObjVal(iter)=h_X;
end
end