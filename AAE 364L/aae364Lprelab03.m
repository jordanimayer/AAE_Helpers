%%%%%
% Jordan Mayer
% AAE 364L
% Prelab 03
% February 13, 2019
%%%%%

clear all; format compact; rehash toolbox;
run('Lab 03 Files/setup_lab_ip01_2_sip.m');
IC_ALPHA0 = 0.2;
disp(A);
disp(B);
options = simset('SrcWorkspace', 'current');

%% Part (i)
close all;

axlim = 18;  % +/- limits for pole plot axes

% pole placement
lambda = [-7+8i -7-8i -5+9i -5-9i]
K = place(A, B, lambda)
sim('Lab 03 Files/s_sip_lqr', [], options);
t = x.time;  % s
x_place = x.signals.values(:,2);  % mm
alpha_place = alpha.signals.values;  % deg

% plot poles
figure();
plot(linspace(-axlim,axlim, 10), zeros(10), '-k');  % real axis
hold on;
plot(zeros(10), linspace(-axlim,axlim,10), '-k');  % imaginary axis
plot(real(lambda), imag(lambda), 'xb');
title('Jordan Mayer - AAE 364L, Prelab 3, Part (i): Poles Used in Pole Placement');
xlabel('Real'); ylabel('Imaginary'); grid on;
axis([-1 1 -1 1] * axlim);

% LQR
Q = diag([0.2 0.2 0.005 0])
R = 0.0001
K = lqr(A, B, Q, R)
sim('Lab 03 Files/s_sip_lqr', [], options);
x_lqr = x.signals.values(:,2);  % mm
alpha_lqr = alpha.signals.values;  % deg

% plot poles
lambda = eig(A - B*K)  % closed-loop poles
figure();
plot(linspace(-axlim,axlim,10), zeros(10), '-k');  % real axis
hold on;
plot(zeros(10), linspace(-axlim,axlim,10), '-k');  % imaginary axis
plot(real(lambda), imag(lambda), 'xb');
title('Jordan Mayer - AAE 364L, Prelab 3, Part (i): Poles Obtained Via LQR');
xlabel('Real'); ylabel('Imaginary'); grid on;
axis([-1 1 -1 1]*axlim);

% plot angle
figure();
plot(t, alpha_place, 'b'); hold on; plot(t, alpha_lqr, 'r');
title('Jordan Mayer - AAE 364L, Prelab 3, Part (i): Angle Results');
xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;
legend('Pole Placement', 'LQR');

% plot position
figure();
plot(t, x_place, 'b'); hold on; plot(t, x_lqr, 'r');
title('Jordan Mayer - AAE 364L, Prelab 3, Part (i): Position Results');
xlabel('Time (s)'); ylabel('Position (mm)'); grid on;
legend('Pole Placement', 'LQR');

% assess results
stats_alpha_place = stepinfo(alpha_place, t, 0);
ts_alpha_place = stats_alpha_place.SettlingTime  % s
stats_x_place = stepinfo(x_place, t, 0);
ts_x_place = stats_x_place.SettlingTime  % s
stats_alpha_lqr = stepinfo(alpha_lqr, t, 0);
ts_alpha_lqr = stats_alpha_lqr.SettlingTime  % s
stats_x_lqr = stepinfo(x_lqr, t, 0);
ts_x_lqr = stats_x_lqr.SettlingTime  % s

%% Part (ii)
close all; clear all;
run('Lab 03 Files/setup_lab_ip01_2_sip.m');
options = simset('SrcWorkspace', 'current');
Ai = [A, zeros(4,1); -1, zeros(1,4)]
Bi = [B; 0]

% obtain gains
Q = diag([0.9 0.9 0 0 2.7])
R = 0.0001
K = lqr(Ai, Bi, Q, R)
sim('Lab 03 Files/aae364pinv2', [], options);
t = x.time;  % s
x_lqr = x.signals.values(:,2);  % mm
alpha_lqr = alpha.signals.values;  % deg

% plot angle
figure();
plot(t, alpha_lqr, 'r');
title('Jordan Mayer - AAE 364L, Prelab 3, Part (ii): Angle Results');
xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;

% plot position
figure();
plot(t, x_lqr, 'r');
title('Jordan Mayer - AAE 364L, Prelab 3, Part (ii): Position Results');
xlabel('Time (s)'); ylabel('Position (mm)'); grid on;