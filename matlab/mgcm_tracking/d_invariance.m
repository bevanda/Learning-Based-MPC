%% Quantiifed uncertainty polytope W
W_lb = [-0.0019; -0.0184; -0.0019; -0.0144; -0.0019; -0.0478; -0.0020; -0.0649];
W_ub = [ 0.0019; 0.0177; 0.0020; 0.0161; 0.0020; 0.0484; 0.0020; 0.0546];
% Form a set with 1000 points taken from the intervals
W=[W_lb(1):(W_ub(1)-W_lb(1))/999:W_ub(1)
W_lb(2):(W_ub(2)-W_lb(2))/999:W_ub(2)
W_lb(3):(W_ub(3)-W_lb(3))/999:W_ub(3)
W_lb(4):(W_ub(4)-W_lb(4))/999:W_ub(4)
W_lb(5):(W_ub(5)-W_lb(5))/999:W_ub(5)
W_lb(6):(W_ub(6)-W_lb(6))/999:W_ub(6)
W_lb(7):(W_ub(7)-W_lb(7))/999:W_ub(7)
W_lb(8):(W_ub(8)-W_lb(8))/999:W_ub(8)];

% Calculate the invariant sets
H(:,:,1) = [min(W'); max(W')];
for ii=2:1:N
C = ((A + (B*K_gain))^(ii-1))*W; C=C';
cc = H(:,:,ii-1);
cc_lb = cc(1,:);
cc_ub = cc(2,:);
CC=[cc_lb(1):(cc_ub(1)-cc_lb(1))/999:cc_ub(1)
cc_lb(2):(cc_ub(2)-cc_lb(2))/999:cc_ub(2)
cc_lb(3):(cc_ub(3)-cc_lb(3))/999:cc_ub(3)
cc_lb(4):(cc_ub(4)-cc_lb(4))/999:cc_ub(4)
cc_lb(5):(cc_ub(5)-cc_lb(5))/999:cc_ub(5)
cc_lb(6):(cc_ub(6)-cc_lb(6))/999:cc_ub(6)
cc_lb(7):(cc_ub(7)-cc_lb(7))/999:cc_ub(7)
cc_lb(8):(cc_ub(8)-cc_lb(8))/999:cc_ub(8)];
CC=CC';
[Msum,MsumRange]=MinkSum(C,CC);
H(:,:,ii) =MsumRange;
cc=[]; C=[]; CC=[];
% Calculate the bounds for X and U
X_set_horizon=[]; U_set_horizon=[];
for ii=1:1:N
cc = H(:,:,ii);
cc_lb = cc(1,:);
cc_ub = cc(2,:);
CC=[cc_lb(1):(cc_ub(1)-cc_lb(1))/999:cc_ub(1)
cc_lb(2):(cc_ub(2)-cc_lb(2))/999:cc_ub(2)
cc_lb(3):(cc_ub(3)-cc_lb(3))/999:cc_ub(3)
cc_lb(4):(cc_ub(4)-cc_lb(4))/999:cc_ub(4)
cc_lb(5):(cc_ub(5)-cc_lb(5))/999:cc_ub(5)
cc_lb(6):(cc_ub(6)-cc_lb(6))/999:cc_ub(6)
cc_lb(7):(cc_ub(7)-cc_lb(7))/999:cc_ub(7)
cc_lb(8):(cc_ub(8)-cc_lb(8))/999:cc_ub(8)];
CC=CC';
PdiffRange=PontDiff(X_set',CC);
X_set_horizon(:,:,ii)=PdiffRange;
PdiffRange=[];
KCC=K_gain*CC';
PdiffRange=PontDiff(U_set',KCC);
U_set_horizon(:,:,ii)=PdiffRange;
end