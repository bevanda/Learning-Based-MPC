close all; addpath('data/casadi/');
%% Plot the solving time througout iteration AMD
load('solve_sample_full.mat');% Create Dependent Variable ‘Experiments’ Data
y1 = [solve_times_1; solve_times_2; solve_times_3; solve_times_4; solve_times_5];   
load('solve_sample_reg.mat');
y2 = [solve_times_1; solve_times_2; solve_times_3; solve_times_4; solve_times_5];  

x = 1:500;                                          % Create Independent Variable                           
N = size(y1,1);                                      % Number of ‘Experiments’ In Data Set
CI95 = tinv([0.025 0.975], N-1); % Calculate 95% Probability Intervals Of t-Distribution

y1Mean = mean(y1);                                    % Mean Of All Experiments At Each Value Of ‘x’
y1SEM = std(y1)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
                   
y1CI95 = bsxfun(@times, y1SEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
figure;
plot(x, y1Mean,'r','LineWidth',2.0)                                      % Plot Mean Of All Experiments
hold on
plot(x, y1CI95+y1Mean, '-.r'); hold on;                              % Plot 95% Confidence Intervals Of All Experiments

y2Mean = mean(y2);                                    % Mean Of All Experiments At Each Value Of ‘x’
y2SEM = std(y2)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
y2CI95 = bsxfun(@times, y2SEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
plot(x, y2Mean,'b','LineWidth',2.0);                                      % Plot Mean Of All Experiments
hold on;
plot(x, y2CI95+y2Mean, '-.b'); hold on;                              % Plot 95% Confidence Intervals Of All Experiments
hold off;
grid;
xlabel('iterations');
ylabel('solve time');
legend('trackingNMPC','CI95','CI95','regularNMPC','CI95','CI95')
title('NMPC solve times for a 500 iteration control window');
%% Plot the solving time througout iteration Intel
load('intelCPU_solve_sample_full.mat');% Create Dependent Variable ‘Experiments’ Data
y3 = [solve_times_1; solve_times_2; solve_times_3; solve_times_4; solve_times_5];   
% load('solve_sample_reg.mat');
% y4 = [solve_times_1; solve_times_2; solve_times_3; solve_times_4; solve_times_5];  

x = 1:500;                                          % Create Independent Variable                           
N = size(y3,1);                                      % Number of ‘Experiments’ In Data Set
CI95 = tinv([0.025 0.975], N-1); % Calculate 95% Probability Intervals Of t-Distribution

y3Mean = mean(y3);                                    % Mean Of All Experiments At Each Value Of ‘x’
y3SEM = std(y3)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
                   
y3CI95 = bsxfun(@times, y3SEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
figure;
plot(x, y3Mean,'r','LineWidth',2.0)                                      % Plot Mean Of All Experiments
hold on
plot(x, y3CI95+y3Mean, '-.r'); hold on;                              % Plot 95% Confidence Intervals Of All Experiments

% y4Mean = mean(y4);                                    % Mean Of All Experiments At Each Value Of ‘x’
% y4SEM = std(y4)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
% y4CI95 = bsxfun(@times, y4SEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
% plot(x, y4Mean,'b','LineWidth',2.0);                                      % Plot Mean Of All Experiments
% hold on;
% plot(x, y4CI95+y4Mean, '-.b'); hold on;                              % Plot 95% Confidence Intervals Of All Experiments
hold off;
grid;
xlabel('iterations');
ylabel('solve time');
legend('trackingNMPC','CI95','CI95')
title({'NMPC solve times for a 500 iteration control window', '(Intel CPU)'});
%% Plot the distribution of the measurements
figure; 
histfit(y3(:)); hold on;  
% histfit(y4(:)); hold on; 
plot(mean(y3(:)), 0,'x','LineWidth', 2.0); hold on;
% plot(mean(y4(:)), 0,'x','LineWidth', 2.0); hold on;
xlabel('solve times');
ylabel('number of samples');
% legend('tNMPC hist','∼N(µ, σ^2 )','rNMPC hist','∼N(µ, σ^2 )','mean tNMPC','mean rNMPC');
legend('tNMPC hist','∼N(µ, σ^2 )','mean tNMPC');
title('Solve time histogram comparable Intel CPU');

%% Plot the solving time througout iteration Intel
load('intelCPU_solve_sample_fullLMPC.mat');% Create Dependent Variable ‘Experiments’ Data
y3 = [solve_times_1; solve_times_2; solve_times_3; solve_times_4; solve_times_5];   
% load('solve_sample_reg.mat');
% y4 = [solve_times_1; solve_times_2; solve_times_3; solve_times_4; solve_times_5];  

x = 1:500;                                          % Create Independent Variable                           
N = size(y3,1);                                      % Number of ‘Experiments’ In Data Set
CI95 = tinv([0.025 0.975], N-1); % Calculate 95% Probability Intervals Of t-Distribution

y3Mean = mean(y3);                                    % Mean Of All Experiments At Each Value Of ‘x’
y3SEM = std(y3)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
                   
y3CI95 = bsxfun(@times, y3SEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
figure;
plot(x, y3Mean,'r','LineWidth',2.0)                                      % Plot Mean Of All Experiments
hold on
plot(x, y3CI95+y3Mean, '-.r'); hold on;                              % Plot 95% Confidence Intervals Of All Experiments

% y4Mean = mean(y4);                                    % Mean Of All Experiments At Each Value Of ‘x’
% y4SEM = std(y4)/sqrt(N);                              % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
% y4CI95 = bsxfun(@times, y4SEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ‘x’
% plot(x, y4Mean,'b','LineWidth',2.0);                                      % Plot Mean Of All Experiments
% hold on;
% plot(x, y4CI95+y4Mean, '-.b'); hold on;                              % Plot 95% Confidence Intervals Of All Experiments
hold off;
grid;
xlabel('iterations');
ylabel('solve time');
legend('trackingLMPC','CI95','CI95')
title({'NMPC solve times for a 500 iteration control window', '(Intel CPU)'});
%% Plot the distribution of the measurements
figure; 
histfit(y3(:)); hold on;  
% histfit(y4(:)); hold on; 
plot(mean(y3(:)), 0,'x','LineWidth', 2.0); hold on;
% plot(mean(y4(:)), 0,'x','LineWidth', 2.0); hold on;
xlabel('solve times');
ylabel('number of samples');
% legend('tNMPC hist','∼N(µ, σ^2 )','rNMPC hist','∼N(µ, σ^2 )','mean tNMPC','mean rNMPC');
legend('tLMPC hist','∼N(µ, σ^2 )','mean tNMPC');
title('Solve time histogram comparable Intel CPU');