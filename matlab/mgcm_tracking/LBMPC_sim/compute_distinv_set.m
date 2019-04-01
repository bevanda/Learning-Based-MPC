function Z_approx = compute_distinv_set(Ak, W, n_order, alpha)
    % W: Polyhedron of system noise
    % We could obtain dist_inv_set Z by computing an infinite geometric series,
    %  which is not practicall to get. So, we approximate this by trancating the polynomial.
    Z_approx = W;
    for n = 1:n_order
        Z_approx = Z_approx + Ak^n*W;
    end
    Z_approx = Z_approx*((1-alpha)^-1);
    % which takes the form of Z = ((1-alpha)^-1)*(W + Ak*W + Ak^2*W + ... Ak^n_ordr*W).
    % where + denotes Minkowski addition.
end