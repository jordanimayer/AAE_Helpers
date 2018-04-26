%%%%%
% Jordan Mayer
% AAE 364
% HW 12
%
% Perform helpful calculations and create Bode and Nyquist plots for
% verification of hand-drawn sketches.
%%%%%

close all;      % close any open figures

% Problem 1: B-7-13
% Bode plot to make hand-drawn Nyquist plot
num = 1;
den = [1 0.2 1 1];
figure('Units', 'inches', 'Position', [1 4 5 4]);
bode(num, den); grid on;
title('Bode Plot for Problem 1: B-7-13');
% Nyquist plot to confirm hand-drawn sketch
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]);
nyquist(G); axis equal;
title('Nyquist Plot for Problem 1: B-7-13');
% Poles of G(s) to assess stability
P = pole(G);
fprintf('\nProblem 1: B-7-13\nP = \n');
disp(P);

% Problem 1: B-7-15
% Bode plot to make hand-drawn Nyquist plot
num = 1;
den = [1 -1 0];
figure('Units', 'inches', 'Position', [1 4 5 4]);
bode(num, den); grid on;
title('Bode Plot for Problem 1: B-7-15');
% Nyquist plot to confirm hand-drawn sketch
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]);
nyquist(G); axis equal;
title('Nyquist Plot for Problem 1: B-7-15');

% Problem 1: B-7-24
% Bode plot to confirm hand-drawn sketch
num = 25;
den = [1 11 10 0];
figure('Units', 'inches', 'Position', [1 4 5 4]);
bode(num, den); grid on;
title('Bode Plot for Problem 1: B-7-24');
% Gain crossover frequency to find phase margin
w_gc = roots([1 0 101 0 100 0 -625]);
fprintf('\nProblem 1:B-7-24\nPossible gain crossover frequencies:\n');
disp(w_gc);

% Problem 2
num = [1.1057 0.1900];
den = [1 0.7385 0.8008 0];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]);
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2');
