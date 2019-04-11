function [xk1, uk] = getTransitionsTrue(xk,xw,r0,K)
%% Closed loop system with prestabilisation /w pole placement
%
% xk1 is the states at time k+1.
uk = -K*(xk-xw)+r0;
[xk1] = trueDynamics(xk, uk);

end
