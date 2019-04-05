%%%%%
% Jordan Mayer
% AAE 364L
% Lab 03
% February 27, 2019
%
% Note that, in experimental results, angle is labeled as theta instead of
% alpha.
%%%%%

close all; clear all; format compact; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 02\Lab 02 Files');
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 02\Lab 02 Experiments');

%% Part (i)
close all;

load('lab2_1_1.mat');
t = part1_theta.time;
alpha = rad2deg(part1_theta.signals.values);

% plot experimental results
figure();
plot(t, alpha, '-b'); xlabel('Time (s)'); ylabel('Angle (deg)');
title('Jordan Mayer - AAE 364L, Lab 03, Part (i): Experimental Results');
grid on;

wp_th = 4.7532;  % theoretical pendulum frequency, rad/s
wp_ex = 4.7842;  % experimental pendulum frequency, rad/s
p_err = abs(wp_th - wp_ex)/wp_th * 100  % percent error, %

%% Part ii
close all; clear all;

load('lab2_2_w7');
t_w3 = part2_w3_theta.time;
t_wn = part2_wn_theta.time;
t_w7 = part2_w7_theta.time;
alpha_w3 = rad2deg(part2_w3_theta.signals.values);
alpha_wn = rad2deg(part2_wn_theta.signals.values);
alpha_w7 = rad2deg(part2_w7_theta.signals.values);

% 3 rad/s results
figure();
plot(t_w3, alpha_w3, '-b'); xlabel('Time (s)'); ylabel('Angle (deg)');
title('Jordan Mayer - AAE 364L, Lab 03, Part (ii): Results for frequency 3 rad/s');
grid on;
mag_w3 = deg2rad(12.22)/3  % transfer function magnitude
mag_w3_dB = 20*log10(mag_w3)  % transfer function magnitude, decibels

% natural frequency results
figure();
plot(t_wn, alpha_wn, '-b'); xlabel('Time (s)'); ylabel('Angle (deg)');
title('Jordan Mayer - AAE 364L, Lab 03, Part (ii): Results for natural frequency');
grid on;
mag_wn = deg2rad(57.22)/3  % transfer function magnitude
mag_wn_dB = 20*log10(mag_wn)  % transfer function magnitude, decibels

% 7 rad/s results
figure();
plot(t_w7, alpha_w7, '-b'); xlabel('Time (s)'); ylabel('Angle (deg)');
title('Jordan Mayer - AAE 364L, Lab 03, Part (ii): Results for frequency 7 rad/s');
grid on;
mag_w7 = deg2rad(13.54)/3  % transfer function magnitude
mag_w7_dB = 20*log10(mag_w7)  % transfer function magnitude, decibels