xl = load('DMS_tLBMPC_q500_pretrained.mat');
xll = load('DMS_tLMPC_K.mat'); load('DMS_tLBMPC_q500.mat');
[l,n]=size(xl.xlo);
figure;
    for i=1:l
        subplot(l,1,i);
        plot(xlo(i,:),'Linewidth',1.5); hold on;
        plot(xll.xl(i,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
        plot(xl.xlo(i,:),'Linewidth',1.5,'LineStyle','--','Color','g'); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['x_{',num2str(i),'}'];
        else
            lableTex = ['u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
        legend('LBMPC online learn','LMPC','LBMPC offline learn')
    end
    
  %%  
xl = load('DMS_tLBMPC_q500_pretrained.mat');
xll = load('tNMPC.mat'); load('DMS_tLBMPC_q500.mat');
[l,n]=size(xl.xlo);
figure;
    for i=1:l
        subplot(l,1,i);
        plot(xlo(i,:),'Linewidth',1.5); hold on;
        plot(xll.xnl(i,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
        plot(xl.xlo(i,:),'Linewidth',1.5,'LineStyle','--','Color','g'); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['x_{',num2str(i),'}'];
        else
            lableTex = ['u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
        legend('LBMPC online learn','NMPC','LBMPC offline learn')
    end
    
    %%
load('DSS_tLMPC.mat'); 
lll=load('DMS_tLMPC_K.mat');
[l,n]=size(xl);
figure;
    for i=1:l
        subplot(l,1,i);
        plot(xl(i,:),'Linewidth',1.5); hold on;
        plot(lll.xl(i,:),'Linewidth',1.5,'LineStyle','--'); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['x_{',num2str(i),'}'];
        else
            lableTex = ['u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
        legend('LMPC','LMPC ')
    end
    %%
load('DMS_tLBMPC.mat');
load('DMS_tLMPC_K.mat'); 
load('tNMPC.mat');
[l,n]=size(xl);
figure;
    for i=1:l
        subplot(l,1,i);
        plot(xlo(i,:),'Linewidth',1.5); hold on;
        plot(xl(i,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
        plot(xnl(i,:),'Linewidth',1.5,'LineStyle','--','Color','g'); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['x_{',num2str(i),'}'];
        else
            lableTex = ['u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
        legend('LBMPC','LMPC','NMPC')
    end