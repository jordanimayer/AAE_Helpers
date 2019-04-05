%%%%%
% Jordan Mayer
% AAE 364L
% Lab 04
% March 20, 2019
%
% Perform some basic arithmetic to obtain results. Create plots using both
% experimental data and simulations.
%%%%%

%% Preliminary setup

close all; clear all; format compact; format long; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 04\Lab 04 Files');
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 04\Lab 04 Data');

%% Part (i)

load('lab4_part2');  % load experimental results
m_s = 0.156;  % kg
m_h = 1.17;  % kg
x = x/100;  % m
x_1 = x1/100;  % m
p_0 = -p0;  % deg
R_c0 = m_s*x_1/m_h  % m
m = m_s + m_h  % kg
R_c = m_s * (x_1 - x)/m  % m
h = R_c/tan(abs(p_0))

%% Part (ii)

load('lab4_part3');  % load experimental results

% extract results
t = theta_average.time;  % s
t_0 = 12.59;  % s, time of helicopter release
after_start_indices = find(t > t_0);
start_index = after_start_indices(1);
t = t(start_index:end) - t_0;  % s
theta = theta_average.signals.values;  % deg
theta = theta(start_index:end);  % deg

% plot pitch angle results
figure();
plot(t, theta, '-b'); grid on;
title('Jordan Mayer - AAE 364L, Lab 04, Part (ii): Pitch Angle Results');
xlabel('Time, s'); ylabel('Pitch Angle, deg');

% calculate omega and sigma
omega = 2*pi*4/13.24  % rad/s
n = 5;  % number of peaks
theta_1 = theta(1)  % deg, first peak
theta_5 = 6.491  % deg, fifth peak
sigma = omega/(2*pi*(n-1)) * log(theta_1/theta_5)  % rad/s

% plot pitch angle results along with enve
env = theta_1 * exp(-sigma*t);  % amplitude envelope, deg
figure();
plot(t, theta, '-b'); hold on; plot(t, env, '--k'); plot(t, -env, '--k');
title('Jordan Mayer - AAE 364L, Lab 04, Part (ii): Pitch Angle Results');
xlabel('Time, s'); ylabel('Pitch Angle, deg'); grid on;
legend('Experimental results', 'Expected envelope');

g = 9.81  % gravitational acceleration, m/s^2
J_p_prime = m*g*h/(omega^2 + sigma^2)  % kg*m^2
J_p = J_p_prime - m_s*(x_1^2 - x^2)  % kg*m^2
c_p = 2*sigma*J_p_prime  % kg/s

% obtain simulated results
J_p_prime = 0.0181;
c_p = 0.00225;
sim('aae364L_lab04_part3');
t_sim = theta_sim.time;  % s
theta_sim = theta_sim.signals.values;  % deg

% plot simulated results
figure();
plot(t, theta, '-b'); hold on; plot(t_sim, theta_sim, '-r');
plot(t, env, '--k'); plot(t, -env, '--k');
title('Jordan Mayer - AAE 364L, Lab 04, Part (ii): Pitch Angle Results');
xlabel('Time, s'); ylabel('Pitch Angle, deg'); grid on;
legend('Experimental results', 'Simulated results', 'Expected envelope');

%% Part (iii), clockwise
load('lab4_part4_cw');  % load experimental data

% extract experimental data
t = psi_dot_average.time;  % s
t_start = 3.274;  % time at which helicopter was pushed, s
t_end = 8.81;  % time at which helicopter stopped, s
ind_after_start = find(t > t_start);
ind_start = ind_after_start(1);
ind_after_end = find(t > t_end);
ind_end = ind_after_end(1);
t = t(ind_start:ind_end) - t_start;  % s
psi_dot = psi_dot_average.signals.values;  % deg
psi_dot = psi_dot(ind_start:ind_end);  % deg

t = t(1:round((0.75*end)));  % only first 75% of data, s
psi_dot = psi_dot(1:round((0.75*end)));  % only first 75% of data, deg
psi_dot = abs(psi_dot);  % must make positive, otherwise logarithms break

% calculate sigma and c_y
Z = log(psi_dot/psi_dot(1));
sigma = -inv(t' * t) * t' * Z
J_shaft = 0.0039;  % kg*m^2
c_y = sigma*(J_p + J_shaft)  % kg*m^2/s

% plot experimental and theoretical results
figure();
plot(t, psi_dot, '-b'); hold on; plot(t, psi_dot(1)*exp(-sigma*t), '-r');
title('Jordan Mayer, AAE 364L, Lab 04, Part (iii), Clockwise Results');
xlabel('Time, s'); ylabel('Yaw rate, deg/s'); grid on;
legend('Experimental results', 'Theoretical results');

%% Part (iii), counterclockwise
load('lab4_part4_ccw');  % load experimental data

% extract experimental data
t = psi_dot_average.time;  % s
t_start = 3.274;  % time at which helicopter was pushed, s
t_end = 8.81;  % time at which helicopter stopped, s
ind_after_start = find(t > t_start);
ind_start = ind_after_start(1);
ind_after_end = find(t > t_end);
ind_end = ind_after_end(1);
t = t(ind_start:ind_end) - t_start;  % s
psi_dot = psi_dot_average.signals.values;  % deg
psi_dot = psi_dot(ind_start:ind_end);  % deg
t = t(1:round((0.75*end)));  % only first 75% of data, s
psi_dot = psi_dot(1:round((0.75*end)));  % only first 75% of data, deg
psi_dot = abs(psi_dot);  % must make positive, otherwise logarithms break

% calculate sigma and c_y
Z = log(psi_dot/psi_dot(1));
sigma = -inv(t' * t) * t' * Z
c_y = sigma*(J_p + J_shaft)  % kg*m^2/s

% plot experimental and theoretical results
figure();
plot(t, psi_dot, '-b'); hold on; plot(t, psi_dot(1)*exp(-sigma*t), '-r');
title('Jordan Mayer, AAE 364L, Lab 04, Part (iii), Counterclockwise Results');
xlabel('Time, s'); ylabel('Yaw rate, deg/s'); grid on;
legend('Experimental results', 'Theoretical results');

%% Part (iv)
load('lab4_part5');  % load experimental data

%% Part (v)
load('lab4_part6');  % load experimental data

%% Part (vi)
load('lab4_part7');  % load experimental data