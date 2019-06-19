load('tLMPC.mat'); load('tNMPC.mat'); load('tLBMPC.mat');
l=size(xl,1);
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
