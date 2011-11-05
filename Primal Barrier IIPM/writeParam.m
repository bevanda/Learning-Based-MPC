% writeParam.m

fid = fopen('ConstrParam.bin','w');

fwrite(fid, kappa_start_PhaseII, 'double');     % double kappa_PhaseII_arg
fwrite(fid, kappa_start_PhaseI, 'double');     % double kappa_PhaseI_arg
fwrite(fid, n_iter_PhaseI, 'int');             % int n_iter_arg
fwrite(fid, n_iter_PhaseII, 'int');             % int n_iter_hat_arg
fwrite(fid, mu, 'double');              % double mu_arg
fwrite(fid, eps_barrier, 'double');     % double eps_barrier_arg
fwrite(fid, eps_nt, 'double');          % double eps_nt_arg
fwrite(fid, eps_normRp, 'double');      % double eps_normRp_arg
fwrite(fid, eps_ls, 'double');          % smallest step size in line search
fwrite(fid, alpha_ls, 'double');        % double alpha_ls_arg
fwrite(fid, beta_ls, 'double');         % double beta_ls_arg
fwrite(fid, reg_PhaseII, 'double');      % regularization term in PhaseI
fwrite(fid, reg_PhaseI, 'double');     % regularization term in PhaseII
fwrite(fid, weight_PhaseI, 'double');    % double weight_hat

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
