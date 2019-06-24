load('DSS_tLMPC.mat'); load('DSS_tNMPC.mat'); load('DMS_tLBMPC_q50.mat');
[l,n]=size(xl);
figure;
    for i=1:l
        subplot(l,1,i);
        plot(xlo(i,:),'Linewidth',1.5); hold on;
        plot(xnl(i,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
        plot(xl(i,:),'Linewidth',1.5,'LineStyle','--'); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['x_{',num2str(i),'}'];
        else
            lableTex = ['u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
        legend('LBMPC','NMPC','LMPC')
    end
