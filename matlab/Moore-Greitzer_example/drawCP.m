
clear;

out = swing_up_RUN();

x1 = out.STATES(:,2);
x2 = out.STATES(:,3);
x3 = out.STATES(:,4);
x4 = out.STATES(:,5);

Ts = 0.1;

Duration = 10;

figure;
subplot(2,2,1);
plot(0:Duration/(length(x1)-1):Duration,x1,'Linewidth',3);
grid on
xlabel('time');
ylabel('z');
title('cart position');
subplot(2,2,2);
plot(0:Duration/(length(x2)-1):Duration,x2,'Linewidth',3);
grid on
xlabel('time');
ylabel('zdot');
title('cart velocity');
subplot(2,2,3);
plot(0:Duration/(length(x3)-1):Duration,x3,'Linewidth',3);
grid on
xlabel('time');
ylabel('theta');
title('pendulum angle');
subplot(2,2,4);
plot(0:Duration/(length(x4)-1):Duration,x4,'Linewidth',3);
grid on
xlabel('time');
ylabel('thetadot');
title('pendulum velocity');