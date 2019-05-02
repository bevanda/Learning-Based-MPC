function [mRPIS] = calc_mRPIS(Ak,F_w,h_w,eps)
% minimal Robust Positively Invariant Set (mRPIS) calculation
%
% Inputs:
%   Ak: transition matrix of the autonomous system
%   F_w,h_w: H-representation of the uncertainty bounds
%   s_0: s-step exact reachability calculation
%   eps: outer \eps-approximation of mRPIS
%
% Output:
%   mRPIS 
%
options = optimset('Display', 'off');
s=0; % init s value
n=size(Ak,1);% system dimension
N=eye(n); % standard basis vectors of R^n
while 1
    s=s+1;
   
    %%%%%%% get \alpha %%%%%%%
    I = length(h_w);
    alphas = zeros(I,1);
    for i = 1:I
        % for each constraint calculate:
        % max_{x}(F_u(i,:)*x) s.t. F_v*x <= h_v
        % support function calculation
        [~, fval] = linprog(-((Ak^s)'*(F_w(i,:))'),F_w,h_w,[], [], [], [], [], options);
        % store alpha values
        alphas(i) = -fval/h_w(i);
    end
    alpha=max(alphas);
    
    %%%%%%% get M %%%%%%%
    M_temp =zeros(n,1);

    for j = 1:n    
        Mp_sum = 0;
        Mm_sum = 0;
        for ind=0:s-2
            % support function calculation
            [~, fvalp] = linprog(-((Ak^ind)'*(N(:,j))),F_w,h_w,[], [], [], [], [], options);
            fval_p=-fvalp;
            [~, fvalm] = linprog(((Ak^ind)'*(N(:,j))),F_w,h_w,[], [], [], [], [], options);
            fval_m=-fvalm;
            
            Mp_sum=Mp_sum+fval_p;
            Mm_sum=Mm_sum+fval_m;
        end
        % store values
        M_temp(j)=max(Mp_sum,Mm_sum);
    end
    M=max(M_temp);
    
    %%%%%%% stopping condition %%%%%%%
    stop_cond=(eps)/(eps+M);
    if (alpha <= stop_cond)
        break;
    end
end

polyH = Polyhedron(F_w,h_w); 
polyH.computeVRep();
polyH.minVRep();
polyV=Polyhedron(polyH.V);
Fs=reach_set(Ak,polyV,s);

scale=1.0/(1.0-alpha);
mRPIS=scale*Fs;
end

