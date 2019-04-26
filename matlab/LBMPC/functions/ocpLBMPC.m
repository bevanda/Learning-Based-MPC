function [sysHistory,art_refHistory,true_refHistory]...
          =ocpLBMPC(x,x_wp,x_wp_init,x_wp_ref,u_wp,...
                    N,Ts,iterations,options,opt_var,data,...
                    A,B,Kstabil,Q,R,P,T,Mtheta,LAMBDA,PSI,m,...
                    F_x,h_x,F_u,h_u,F_w_N,h_w_N,F_x_d,h_x_d,...
                    sysHistory,art_refHistory,true_refHistory)
%==========================================================================
% Solving the Learnin-Based Model Predictive optimal control problem 
% at every timestep
%==========================================================================
for k = 1:(iterations)      
    fprintf('iteration no. %d/%d \n',k,iterations);
    if k>1     
        % data acquisition 
        X=[x(1:2)-x_wp(1:2); u-u_wp]; %[δx1;δx2;δu]
        Y=((x_k1-x_wp)-(A*(x-x_wp)+B*(u-u_wp))); %[δx_true-δx_nominal]
       
        x=x_k1; % update state vars for estimation
        q=100; % moving window of q datapoints 
        data=update_data(X,Y,q,k,data); % update data
        
    % get the real state w.r.t. equilibrium
        dx=x-x_wp;
    else
        dx=x_wp_init;
    end
    
    % Solve the OCP
    COSTFUN = @(var) costLBMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        dx,x_wp_ref,N,reshape(var(1:m),m,1),Q,R,P,T,Kstabil,x_wp,u_wp,LAMBDA,PSI,data,Ts);
    CONSFUN = @(var) constraintsLBMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        dx,N,Kstabil,F_x,h_x,F_u,h_u,F_w_N,h_w_N,F_x_d,h_x_d);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;

    % Implement first optimal control move and update plant states.
    [x_k1, u] = transitionTrue(x,c,x_wp,u_wp,Kstabil,Ts); % plant   
    
    % Save state data for plotting w.r.t. work point x_wp
    [x, u]=wp_shift(x,x_wp,u,u_wp);
    his = [x; u]; 
    
    % Save plant states for display.
    sysHistory = [sysHistory his]; %#ok<*AGROW>
    art_refHistory = [art_refHistory art_ref(1:m)];
    true_refHistory = [true_refHistory x_wp_ref];
    
end