addpath('data/casadi/');

%% N=100
load('DMS_tLBMPC_q100.mat');
load('tLMPC.mat'); 
load('tNMPC.mat');
[l,n]=size(xlo);
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
%%  N=50
load('DMS_N50_tLBMPC_q10.mat');
load('DMS_N50_tLMPC.mat'); 
[l,n]=size(xlo);
figure;
    for i=1:l
        subplot(l,1,i);
        plot(xlo(i,:),'Linewidth',1.5); hold on;
        plot(xl(i,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['x_{',num2str(i),'}'];
        else
            lableTex = ['u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
        legend('LBMPC','LMPC')
    end
%%
load('DMS_N50_tLBMPC_q100.mat');
load('DMS_N50_tLMPC.mat'); 
[l,n]=size(xlo);

figure;
    for i=1:l
        subplot(l,1,i);
        plot(xlo(i,:),'Linewidth',1.5); hold on;
        plot(xl(i,:),'Linewidth',1.5,'LineStyle','-.'); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['x_{',num2str(i),'}'];
        else
            lableTex = ['u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
        legend('LBMPC','LMPC')
    end