%==========================================================================
% Write to quad.dat in a binary format.  Note that the transposes/order of
% the variables is important to compensate for the row-major/column-major
% differences in array storage between Matlab, C++, and Fortran
%==========================================================================
function write_quad_dat(fname_dat, fname_mat, ...
    constraint_count, ...
    A, B, K, ...
    Proj_X_t, Proj_U_t, ...
    P, Q, dlqr_controlweight, R, ...
    Aineq, bineq, b_crx, ...
    sub_block, uncertainty_structure, ...
    uncertainty_lower_bounds, uncertainty_upper_bounds, ...
    ALPHA, MU_VEC, d_0)

disp(['Writing to ' fname_dat '...']);

fid = fopen(fname_dat, 'w');

fwrite(fid, constraint_count, 'int32');

fwrite(fid, [A'; B']', 'double'); % becomes M in C++
fwrite(fid, K', 'double');

fwrite(fid, Proj_X_t', 'double');
fwrite(fid, Proj_U_t', 'double');

fwrite(fid, sqrtm(P)', 'double');
fwrite(fid, sqrtm(Q)', 'double');
fwrite(fid, sqrtm(dlqr_controlweight*R)', 'double');

fwrite(fid, Aineq, 'double');
fwrite(fid, bineq, 'double');

fwrite(fid, b_crx', 'double');

fwrite(fid, sub_block, 'double');
fwrite(fid, uncertainty_structure', 'double');

fwrite(fid, uncertainty_lower_bounds', 'double');
fwrite(fid, uncertainty_upper_bounds', 'double');

fwrite(fid, ALPHA, 'double');
fwrite(fid, MU_VEC, 'double');

fwrite(fid, d_0, 'double');

fclose(fid);

save(fname_mat, 'A', 'B', 'd_0')




