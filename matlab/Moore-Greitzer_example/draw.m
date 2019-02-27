%% Run
out = mgcm_RUN();

u = out.STATES(:,1);
x1 = out.STATES(:,2);
x2 = out.STATES(:,3);

Ts = 0.2;

Duration = 10;

figure;
subplot(2,2,1);
plot(0:Duration/(length(x1)-1):Duration,x1,'Linewidth',3);
grid on
xlabel('time');
ylabel('x1');
title('cart position');
subplot(2,2,2);
plot(0:Duration/(length(x2)-1):Duration,x2,'Linewidth',3);
grid on
xlabel('time');
ylabel('x1');
title('cart velocity');
subplot(2,2,3);
plot(0:Duration/(length(u)-1):Duration,u,'Linewidth',3);
grid on
xlabel('time');
ylabel('theta');
title('Control');
subplot(2,2,4);
plot(0:x1,x2,'Linewidth',3);
grid on
xlabel('x2');
ylabel('x1');
title('State-Space');