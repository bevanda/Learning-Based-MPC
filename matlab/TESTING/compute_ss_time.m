function t_ss = compute_ss_time(t,x,target, eps)
    i_ss = x > target - eps & x < target + eps;
    i_ss_flip = ~fliplr(i_ss');
    i_last = length(x) - find(i_ss_flip,1) + 1;
    t_ss = t(min(length(t),i_last + 1));

end %function