function control_error = compute_ctrl_error(t,x,x_target, t_ss)
    dt = t(2)-t(1);  
    x = x(t <= t_ss);
    control_error = sum((x - x_target).^2) * dt;
end %function