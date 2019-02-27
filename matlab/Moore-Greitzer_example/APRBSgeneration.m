%APRBS 
Ts = 0.01 %sampling time
training_length = 100;
t_min=0.1; t_max=2.5; %min & max timestep of APRBS signal in seconds
a_min = -5; a_max = 5; %min & max amplitude of Trainings data
tp = tp_input(training_length/Ts,[t_min,t_max],[a_min, a_max],Ts);