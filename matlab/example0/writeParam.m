% writeParam.m

fid = fopen(fileName,'w');

fwrite(fid, n_iter, 'int');             % int n_iter_arg
fwrite(fid, reg, 'double');      % regularization term in PhaseI
fwrite(fid, eps_primal, 'double');      % residuum norm on r_H
fwrite(fid, eps_dual, 'double');      % residuum norm on r_C
fwrite(fid, eps_mu, 'double');          % residuum norm on gap

fwrite(fid, A, 'double');               % Matrix<Type, _n, _n> &A_arg
fwrite(fid, B, 'double');               % Matrix<Type, _n, _m> &B_arg
fwrite(fid, s, 'double');               % Matrix<Type, _n, 1> &s

fwrite(fid, Q_tilde, 'double');         % Matrix<Type, _n, _n> &Q_tilde_arg
fwrite(fid, Q_tilde_f, 'double');       % Matrix<Type, _n, _n> &Q_tilde_f_arg
fwrite(fid, R, 'double');               % Matrix<Type, _m, _m> &R_arg

for i = 1 : N
    fwrite(fid, Fx{i}, 'double');             % Matrix<Type, _nSt, _n> Fx_arg[]
end

for i = 1 : N
    fwrite(fid, fx{i}, 'double');             % Matrix<Type, _nSt, _n> Fx_arg[]
end

for i = 1 : N
    fwrite(fid, Fu{i}, 'double');             % Matrix<Type, Dynamic, _m, 0, _nInp, _m> Fu_arg[]
end

for i = 1 : N
    fwrite(fid, fu{i}, 'double');             % Matrix<Type, Dynamic, 1, 0, _nInp, 1> fu_arg[]
end

fwrite(fid, F_xTheta, 'double');        % Matrix<Type, _nF_xTheta, _n> &F_xTheta_arg
fwrite(fid, F_theta, 'double');         % Matrix<Type, _nF_xTheta, _m> &F_theta_arg
fwrite(fid, f_xTheta, 'double');        % Matrix<Type, _nF_xTheta, 1> &f_xTheta_arg
fwrite(fid, K, 'double');               % Matrix<Type, _m, _n> &K_arg 

fclose(fid);
