function control_energy = compute_ctrl_energy(t,u,t_ss)
    dt = t(2)-t(1);
    u = u(t <= t_ss);
    control_energy = dt * sum(u.^2);
end %function