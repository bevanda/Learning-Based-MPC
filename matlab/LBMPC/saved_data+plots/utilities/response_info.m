function info=response_info(y,u,t,ys)
    T=t;
    plot(t,y); hold on;
    info=stepinfo(y,T,ys,'SettlingTimeThreshold',0.05);
    info.CumError=compute_ctrl_error(y,T,ys,info.SettlingTime);
    info.CumEnergy=compute_ctrl_energy(u,t,info.SettlingTime);
    info.Overshoot=compute_overshoot(y,ys);
end