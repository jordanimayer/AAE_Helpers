%%%%%
% Jordan Mayer
% AAE 364L
% Experiment 1: Cart on Track
% January 30, 2019
%
% Run simulations and create plots of both theoretical and experimental
% data. Generate basic numerical characteristics of closed loop response to
% different controllers.
%%%%%

clear all; close all; options = simset('SrcWorkspace', 'current');

% Part i

load('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 01\Lab 01 Experiments\1_1_no_EMF.mat')

figure();
plot(y.signals.values, 'b'); hold on;  % plot experimental results
title('Jordan Mayer - AAE 364L, Lab 01, Part (i): Cart Results');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;

t = 0:(7610-2829);  % time, ms
B_eq = 6.2186;  % damping due to viscous friction, kg/s
B_eq = B_eq/1000;  % kg/ms
y_dot0 = 0.00078;  % initial velocity, m/ms
y_inf = 0.1346;  % final position, m
m = 1.0731;  % total cart mass, kg
y_th = (y_dot0*m/B_eq) * (1 - exp(-(B_eq*t/m)));
  % theoretical position, m

figure();
plot(y.signals.values(2829:7610) - 0.0087, 'b'); hold on;  % plot experimental results
plot(t, y_th, 'r'); % plot theoretical results
title('Jordan Mayer - AAE 364L, Lab 01, Part (i): Cart Results');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;
legend('Experimental', 'Theoretical');

%% Part ii

clear all; options = simset('SrcWorkspace', 'current');
load('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 01\Lab 01 Experiments\2_1.mat')

r_0 = 0.4;  % m
k = 50;  % V/m
gamma = 1.72;  % s*A/m
c = 13.943;  % kg/s
m = 1.0731;  % kg
sim('aae364L_lab01_part2_noSat', [], options);
sim('aae364L_lab01_part2_withSat', [], options);

figure();
plot(v.signals.values, 'b'); hold on;  % plot experimental results
plot(v_noSat.signals.values, 'r');
  % plot simulation (no saturation) results
plot(v_withSat.signals.values, 'g');
  % plot simulation (with saturation) results
title('Jordan Mayer - AAE 364L, Lab 01, Part (ii), Voltage Results');
xlabel('Time (ms)'); ylabel('Voltage (V)'); grid on;
legend('Experimental', 'Simulation (no saturation)',...
       'Simulation (with saturation)');

figure();
plot(y.signals.values, 'b'); hold on;  % plot experimental results
plot(y_noSat.signals.values, 'r');
  % plot simulation (no saturation) results
plot(y_withSat.signals.values, 'g');
  % plot simulation (with saturation) results
title('Jordan Mayer - AAE 364L, Lab 01, Part (ii), Position Results');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;
legend('Experimental', 'Simulation (no saturation)',...
       'Simulation (with saturation)');

c = 16;  % kg/s
sim('aae364L_lab01_part2_noSat', [], options);
sim('aae364L_lab01_part2_withSat', [], options);
figure();
plot(y.signals.values, 'b'); hold on;  % plot experimental results
plot(y_noSat.signals.values, 'r');
  % plot simulation (no saturation) results
plot(y_withSat.signals.values, 'g');
  % plot simulation (with saturation) results
title('Jordan Mayer - AAE 364L, Lab 01, Part (ii), Position Results (c = 16 kg/s)');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;
legend('Experimental', 'Simulation (no saturation)',...
       'Simulation (with saturation)');
   
c = 10;  % kg/s
sim('aae364L_lab01_part2_noSat', [], options);
sim('aae364L_lab01_part2_withSat', [], options);
figure();
plot(y.signals.values, 'b'); hold on;  % plot experimental results
plot(y_noSat.signals.values, 'r');
  % plot simulation (no saturation) results
plot(y_withSat.signals.values, 'g');
  % plot simulation (with saturation) results
title('Jordan Mayer - AAE 364L, Lab 01, Part (ii), Position Results (c = 10 kg/s)');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;
legend('Experimental', 'Simulation (no saturation)',...
       'Simulation (with saturation)');

%% Part iii

clear all; options = simset('SrcWorkspace', 'current');
load('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 01\Lab 01 Experiments\3_1_1.mat')

figure();
plot(y.signals.values(1:2000));  % plot experimental results
title('Jordan Mayer - AAE 364L, Lab 01, Part (iii), Position Results');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;

clear all; options = simset('SrcWorkspace', 'current');
load('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 01\Lab 01 Experiments\3_1_2.mat')

k = 10;  % V/m
r = 0.4;  % m
gamma = 1.72;  % s*A/m
c = 13.943;  % kg/s
m = 1.0731;  % kg
fc = 0.2683;  % N
sim('aae364L_lab01_part3', [], options);

plot(e.signals.values(1:3000), 'b');  % plot experimental results
hold on;
plot(e_sim.signals.values, 'r');  % plot simulation results
title('Jordan Mayer - AAE 364L, Lab 01, Part (iii), Error Results');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;
legend('Experimental', 'Theoretical');

%% Part v

clear all; options = simset('SrcWorkspace', 'current');

kp = 100;
kd = 5;
ki = 0;
gamma = 1.72;  % s*A/m
m = 1.0731;  % kg
c = 13.943;  % kg/s
sim('aae364L_lab01_part5', [], options);

figure();
plot(y_sim.signals.values, 'b');  % plot simulated results
title('Jordan Mayer - AAE 364L, Lab 01, Part (v), Simulated Results');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;

y = y_sim.signals.values;
t = y_sim.time;
S = stepinfo(y, t)
e = y(size(y,1)) - 0.5

clear all; options = simset('SrcWorkspace', 'current');

load('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 01\Lab 01 Experiments\5_1_Jordan_40_5_5.mat')

figure();
plot(y.signals.values, 'b');  % plot experimental results
title('Jordan Mayer - AAE 364L, Lab 01, Part (v), Results with Experimental Gains');
xlabel('Time (ms)'); ylabel('Position (m)'); grid on;

y_t = y.signals.values;
t = y.time;
S = stepinfo(y_t, t)
e = y_t(size(y_t,1)) - 0.5

clear all; options = simset('SrcWorkspace', 'current');

kp = 40;
kd = 5;
ki = 5;
gamma = 1.72;  % s*A/m
m = 1.0731;  % kg
c = 13.943;  % kg/s
sim('aae364L_lab01_part5', [], options);
hold on;
plot(y_sim.signals.values, 'r');  % plot simulated results
legend('Experimental', 'Theoretical');

y = y_sim.signals.values;
t = y_sim.time;
S = stepinfo(y, t)
e = y(size(y,1)) - 0.5