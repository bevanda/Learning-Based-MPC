function info=response_info(y,u,t,ys)
%     T=repmat(t,size(y,1),1);
    t=repmat(t,size(u,1),1);
    T=t;
    info=stepinfo(y,T,ys);
    info.CumError=compute_ctrl_error(y,T,ys,info.SettlingTime);
    info.CumEnergy=compute_ctrl_energy(u,t,info.SettlingTime);
end