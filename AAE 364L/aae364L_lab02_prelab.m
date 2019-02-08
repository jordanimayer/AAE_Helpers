%%%%%
% Jordan Mayer
% AAE 364L
% Lab 02
% Pre-Lab
% January 30, 2019
%%%%%

clear all; close all; format compact; rehash toolbox;

%% Part (ii)

A = [0 0 1 0; 0 0 0 1; 0 1.5216 -11.6513 0.0049; 0 -26.1093 26.8458 -0.0841];
B = [0 0 1.5304 -3.5261].';
C = [0 1 0 0];  % for output = angle (alpha)
D = 0;
sys_alpha = ss(A,B,C,D);

bode(sys_alpha); grid on;
title('Jordan Mayer - AAE 364L, Prelab 02, Part (ii): Bode Plot of G_{alpha}');


cd 'Lab 02 Files';
%% Part (iii)
fprintf('\nPart (iii):\n\n');

run('setup_lab_ip02_spg.m');
options = simset('SrcWorkspace', 'current');
X0 = [0 0 0 pi/2].';
K = [0 0 0 0];
sim('s_spg_pp_noK', [], options);
K = place(A, B, [-1.8182+1.9067i -1.8182-1.9067i -27 -28])
sim('s_spg_pp_withK', [], options);

figure();
plot(alpha_noK.time, alpha_noK.signals.values, 'b'); hold on;
plot(alpha_withK.time, alpha_withK.signals.values, 'r');
title('Jordan Mayer - AAE 364L, Prelab 02, Part (iii): Angle Results');
xlabel('Time (s)'); ylabel('Angle (deg)'); grid on; xlim([0 15]);
legend('No Control', 'With Control');

stats_withK = stepinfo(alpha_withK.signals.values, alpha_withK.time, 0);
t_s_withK = stats_withK.SettlingTime

%% Part iv
fprintf('\nPart (iv):\n\n');

clear all;
run('setup_lab_ip02_spg.m');
options = simset('SrcWorkspace', 'current');
Ai = [A zeros(4,1); -1 0 0 0 0];
Bi = [B; 0];
K = place(Ai, Bi, [-3+3i -3-3i -8+4i -8-4i -4.5])
sim('aae364gantry2', [], options);

figure();
plot(x_c.time, x_c.signals.values);
title('Jordan Mayer - AAE 364L, Prelab 02, Part (iv): Position Results');
xlabel('Time (s)'); ylabel('Cart Position (m)'); grid on; xlim([0 15]);

stats_x_c = stepinfo(x_c.signals.values, x_c.time, 0.5);
PO_x_c = stats_x_c.Overshoot
t_s_x_c = stats_x_c.SettlingTime

figure();
plot(alpha.time, alpha.signals.values);
title('Jordan Mayer - AAE 364L, Prelab 02, Part (iv): Angle Results');
xlabel('Time (s)'); ylabel('Pendulum Angle (deg)'); grid on; xlim([0 15]);

stats_alpha = stepinfo(alpha.signals.values, alpha.time, 0);
t_s_alpha = stats_alpha.SettlingTime