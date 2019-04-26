function plot2DSS(x1,x2)
    figure;
    plot(x1,x2,'Linewidth',1.5,'LineStyle','-'); hold on;
    grid on;
    xlabel('x_1');
    ylabel('x_2');
    title('State space');
end