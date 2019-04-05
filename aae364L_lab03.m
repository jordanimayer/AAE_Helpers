%%%%%
% Jordan Mayer
% AAE 364L
% Lab 03
% February 27, 2019
%%%%%

format compact; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 03\Lab 03 Files');
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 03\Lab 03 Experiments');
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 03\Lab 03 Simulations');

%% Part (i)

close all; clear all; 

% load results
load('lab3_part1_sim');  % load simulated results
x_place_sim = x_place_sim/1000;  % m (was mm)
x_lqr_sim = x_lqr_sim/1000;  % m (was mm)
load('lab3_part1');  % load experimental results

% extract experimental results with pole placement controller
t_start_place = 4.6;  % approximate time when controller started, s
t_place = theta_part1_jordan.time;  % s
n_after = find(t_place > t_start_place);
n_start = n_after(1);  % index of controller start time
len = size(t_place,1);  % total length of data
t_place = t_place(n_start:len) - t_start_place;
alpha_place = rad2deg(theta_part1_jordan.signals.values);  % rad
alpha_place = alpha_place(n_start:len);
x_place = x_part1_jordan.signals.values;  % m
x_place = x_place(n_start:len);

% extract experimental results with LQR controller
t_start_lqr = 4.2;  % approximate time when controller started, s
t_lqr = theta_part1_lqr_jordan.time;  % s
n_after = find(t_lqr > t_start_lqr);
n_start = n_after(1);  % index of controller start time
len = size(t_lqr,1);  % total length of data
t_lqr = t_lqr(n_start:len) - t_start_lqr;
alpha_lqr = rad2deg(theta_part1_lqr_jordan.signals.values);  % rad
alpha_lqr = alpha_lqr(n_start:len);
x_lqr = x_part1_lqr_jordan.signals.values;  % m
x_lqr = x_lqr(n_start:len);

% plot position
figure();
plot(t_sim, x_place_sim, '-b'); hold on;
plot(t_sim, x_lqr_sim, '-r');
plot(t_place, x_place, '-g');
plot(t_lqr, x_lqr, '-k');
title('Jordan Mayer - AAE 364L, Lab 03, Part (i): Position Results');
xlabel('Time (s)'); ylabel('Cart Position (m)'); grid on;
legend('Pole Placement (Simulated)', 'LQR (Simulated)', ...
       'Pole Placement (Experimental)', 'LQR (Experimental)');

% plot angle
figure();
plot(t_sim, alpha_place_sim, '-b'); hold on;
plot(t_sim, alpha_lqr_sim, '-r');
plot(t_place, alpha_place, '-g');
plot(t_lqr, alpha_lqr, '-k');
title('Jordan Mayer - AAE 364L, Lab 03, Part (i): Angle Results');
xlabel('Time (s)'); ylabel('Pendulum Angle (deg)'); grid on;
legend('Pole Placement (Simulated)', 'LQR (Simulated)', ...
       'Pole Placement (Experimental)', 'LQR (Experimental)');
 
%% Part (ii)

% load results
load('lab3_part2_sim');  % load simulated results
x_lqr_sim = x_lqr_sim/1000;  % m (was mm)
load('lab3_part2');  % load experimental results

% extract experimental results
t_start_lqr = 8.1;  % approximate time when controller started, s
t_lqr = theta_part3_Jordan.time;  % s
n_after = find(t_lqr > t_start_lqr);
n_start = n_after(1);  % index of controller start time
len = size(t_lqr,1);  % total length of data
t_lqr = t_lqr(n_start:len) - t_start_lqr;
alpha_lqr = rad2deg(theta_part3_Jordan.signals.values);  % rad
alpha_lqr = alpha_lqr(n_start:len);
x_lqr = x_part3_Jordan.signals.values;  % m
x_lqr = x_lqr(n_start:len);

% plot position
figure();
plot(t_sim, x_lqr_sim, '-b'); hold on; plot(t_lqr, x_lqr, '-r');
title('Jordan Mayer - AAE 364L, Lab 03, Part (ii): Position Results');
xlabel('Time (s)'); ylabel('Cart Position (m)'); grid on;
legend('Simulated', 'Experimental');

% plot angle
figure();
plot(t_sim, alpha_lqr_sim, '-b'); hold on; plot(t_lqr, alpha_lqr, '-r');
title('Jordan Mayer - AAE 364L, Lab 03, Part (ii): Angle Results');
xlabel('Time (s)'); ylabel('Pendulum Angle (deg)'); grid on;
legend('Simulated', 'Experimental');