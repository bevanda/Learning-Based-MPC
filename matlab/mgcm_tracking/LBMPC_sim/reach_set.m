function Z_approx = reach_set(Ak, W, n_order)
    % W: Polyhedron of system noise
    % We could obtain dist_inv_set Z by computing an infinite geometric series,
    %  which is not practicall to get. So, we approximate this by trancating the polynomial.
    W.computeVRep();
    Z_approx = W;
    for n = 1:n_order
        Z_approx = Z_approx + Ak*W;
    end
    % which takes the form of Z = (W + Ak*W + Ak^2*W + ... Ak^n_ordr*W).
    % where + denotes Minkowski addition.
end