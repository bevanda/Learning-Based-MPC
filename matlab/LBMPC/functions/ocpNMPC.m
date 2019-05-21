function [sysHistory,art_refHistory,true_refHistory]...
          =ocpNMPC(x,x_wp,x_wp_ref,u_wp,...
                    N,Ts,iterations,options,opt_var,...
                    Kstabil,Q,R,P,T,Mtheta,LAMBDA,PSI,m,...
                    F_x,h_x,F_u,h_u,F_w_N,h_w_N,...
                    sysHistory,art_refHistory,true_refHistory)

% Solving the Nonlinear Model Predictive optimal control problem 
% at every timestep

for k = 1:(iterations)      
    fprintf('iteration no. %d/%d \n',k,iterations);
    % get the true state

    % Solve the OCP   
    COSTFUN = @(var) costNMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        x,x_wp_ref,x_wp,u_wp,N,Ts,reshape(var(1:m),m,1),....
        Q,R,P,T,Kstabil,LAMBDA,PSI);
    CONSFUN = @(var) constraintsNMPC(reshape(var(1:end-m),m,N),reshape(var(end-m+1:end),m,1),...
        x,N,Ts,x_wp,u_wp,Kstabil,...
        F_x,h_x,F_u,h_u,F_w_N,h_w_N);
    opt_var = fmincon(COSTFUN,opt_var,[],[],[],[],[],[],CONSFUN,options);    
    theta_opt = reshape(opt_var(end-m+1:end),m,1);
    c = reshape(opt_var(1:m),m,1);
    art_ref = Mtheta*theta_opt;

    % Implement first optimal control move and update plant states.
    [x, u] = transitionTrue(x,c,x_wp,u_wp,Kstabil,Ts); % plant   
    
    % Save state data for plotting w.r.t. work point x_wp
    his = [x-x_wp; u-u_wp]; 
    
    % Save plant states for display.
    sysHistory = [sysHistory his]; %#ok<*AGROW>
    art_refHistory = [art_refHistory art_ref(1:m)];
    true_refHistory = [true_refHistory x_wp_ref];
    
end