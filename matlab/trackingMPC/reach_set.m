function Z = reach_set(Ak, W, s_order)
    % W: Polyhedron of system noise
    % We could obtain dist_inv_set Z by computing an infinite geometric series,
    %  which is not practicall to get. So, we approximate this by trancating the polynomial.
%     W.computeVRep();
    Z = W;
    for n = 1:s_order-1
        Z = Z + Ak^n*W;
        Z.minVRep();
    end
    % which takes the form of Z = (W + Ak*W + Ak^2*W + ... Ak^n_ordr*W).
    % where + denotes Minkowski addition.
end