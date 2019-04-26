function [x_out, u_out]=wp_shift(x,x_w,u,u_w)
% Shifts (x,u) w.r.t. to working/linearization point (x_w,u_w)
    x_out=x-x_w;
    u_out=u-u_w;
end