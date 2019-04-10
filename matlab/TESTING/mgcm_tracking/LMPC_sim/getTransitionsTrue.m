function [xk1, uk] = getTransitionsTrue(xk,ck,xw,r0,K,Ts)
%% Closed loop system with prestabilisation /w pole placement
%
% xk1 is the states at time k+1.
uk = K*(xk-xw)+ck+r0;
[xk1] = trueDynamics(xk, uk,Ts);

end
