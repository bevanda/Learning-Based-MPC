close all;
%% Plot the solveing time througout iteration
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
%% Plot the distibution of the measurements
figure; 
histfit(y1(:)); hold on;  
histfit(y2(:)); hold on; 
plot(mean(y1(:)), 0,'x','LineWidth', 2.0); hold on;
plot(mean(y2(:)), 0,'x','LineWidth', 2.0); hold on;
xlabel('solve times');
ylabel('number of samples');
legend('tNMPC hist','∼N(µ, σ^2 )','rNMPC hist','∼N(µ, σ^2 )','mean tNMPC','mean rNMPC');
title('Solve time histogram');