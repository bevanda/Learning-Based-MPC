function plotRESPONSE(sysHistory, art_refHistory, t_vec, n, m)
figure;
    for i=1:n+m
        subplot(n+m,1,i);
        plot(t_vec,sysHistory(i,:),'Linewidth',1.5); hold on;
        grid on
        xlabel('time [s]');
        if i<=n
            lableTex = ['\delta x_{',num2str(i),'}'];
        else
            lableTex = ['\delta u_{',num2str(i-n),'}'];
        end
        ylabel(lableTex,'Interpreter','tex' );
    end
end