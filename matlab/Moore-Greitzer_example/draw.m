%% Run
out = mgcm_RUN();

u = out.STATES(:,1);
x1 = out.STATES(:,2);
x2 = out.STATES(:,3);
x3 = out.STATES(:, 4);
x4 = out.STATES(:, 5);

Ts = 0.01;

Duration = 600;

figure;
subplot(2,2,1);
plot(0:Duration/(length(x1)-1):Duration,x1,'Linewidth',3);
grid on
xlabel('time');
ylabel('x1');
title('mass flow');
subplot(2,2,2);
plot(0:Duration/(length(x2)-1):Duration,x2,'Linewidth',3);
grid on
xlabel('time');
ylabel('x1');
title('pressure rise');
subplot(2,2,3);
plot(0:Duration/(length(u)-1):Duration,u,'Linewidth',3);
grid on
xlabel('time');
ylabel('U input');
title('Control');
subplot(2,2,4);
plot(x2,x1,'Linewidth',3);
grid on
xlabel('x1');
ylabel('x2');
title('State-Space');