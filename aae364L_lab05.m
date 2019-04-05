%%%%%
% Jordan Mayer
% AAE 364L
% Lab 05
% April 5, 2019
%%%%%

%% Preliminary setup
close all; clear all; format compact; format short; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 05\Lab 05 Files');
addpath('W:\Courses spring 2019\AAE 364L\AAE 364L Labs\Lab 05\Lab 05 Data');

%% Part (i), Pitch
fprintf('\nPart (i):\n');

fprintf('\nPitch performance:\n\n');
load('part7_15_mayer15');  % load pitch control data
t = (0:9900)/100;  % time, s

% create plots
figure();
plot(t, pitch(:,1), '-b'); hold on; plot(t, pitch(:,2), '-r'); grid on;
title('Jordan Mayer - AAE 364L, Lab 05, Part (i): Pitch Control Results');
xlabel('Time, sec'); ylabel('Pitch Angle, deg');
legend('Desired', 'Actual'); axis([0 99 -35 35]);

figure();
plot(t, yaw(:,1), '-b'); hold on; plot(t, yaw(:,2), '-r'); grid on;
title('Jordan Mayer - AAE 364L, Lab 05, Part (i): Pitch Control Results');
xlabel('Time, sec'); ylabel('Yaw Angle, deg');
legend('Desired', 'Actual'); axis([0 99 -5 5]);

figure();
plot(t, voltage(:,1), '-b'); hold on; plot(t, voltage(:,2), '-r'); grid on;
title('Jordan Mayer - AAE 364L, Lab 05, Part (i): Pitch Control Results');
xlabel('Time, sec'); ylabel('Voltage, V');
legend('Pitch', 'Yaw');

% get pitch performance info
info = stepinfo(pitch(:,2), t);
t_r = info.RiseTime  % s
over = info.Overshoot % deg
t_s = info.SettlingTime  % s
ind_s = find(abs(t - t_s) < 0.005);  % index of settling
e_ss = pitch(ind_s) - 15*sin(2*pi*(1/30)*t_s)  % steady-state error, deg

%% Part (i), Yaw
fprintf('\nYaw performance:\n\n');

load('part7_30_mayer15');  % load pitch control data
t = (0:9900)/100;  % time, s

% create plots
figure();
plot(t, pitch(:,1), '-b'); hold on; plot(t, pitch(:,2), '-r'); grid on;
title('Jordan Mayer - AAE 364L, Lab 05, Part (i): Yaw Control Results');
xlabel('Time, sec'); ylabel('Pitch Angle, deg');
legend('Desired', 'Actual'); axis([0 99 -10 10]);

figure();
plot(t, yaw(:,1), '-b'); hold on; plot(t, yaw(:,2), '-r'); grid on;
title('Jordan Mayer - AAE 364L, Lab 05, Part (i): Yaw Control Results');
xlabel('Time, sec'); ylabel('Yaw Angle, deg');
legend('Desired', 'Actual'); axis([0 99 -40 40]);

figure();
plot(t, voltage(:,1), '-b'); hold on; plot(t, voltage(:,2), '-r'); grid on;
title('Jordan Mayer - AAE 364L, Lab 05, Part (i): Yaw Control Results');
xlabel('Time, sec'); ylabel('Voltage, V');
legend('Pitch', 'Yaw');

% get yaw performance info
info = stepinfo(yaw(:,2), t);
t_r = info.RiseTime  % s
over = info.Overshoot % deg
t_s = info.SettlingTime  % s
ind_s = find(abs(t - t_s) < 0.005);  % index of settling
e_ss = pitch(ind_s) - 30*sin(2*pi*.02*t_s)  % steady-state error, deg