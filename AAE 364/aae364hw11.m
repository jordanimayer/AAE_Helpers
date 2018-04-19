%%%%%
% Jordan Mayer
% AAE 364
% HW 11
%
% Create Bode plots and Nyquist plots to check hand-drawn sketches.
%%%%%

format compact; % use compact formatting for outputs
close all;      % close any open figures

% Problem 1.1
fprintf('\nProblem 1.1\n');
num = [1.1057 0.1900];
den = [1 0.7385 0.8008 0];
G = tf(num, den);
z = zero(G);
p = pole(G);
fprintf('Zeros of G(s):\n');
disp(z);
fprintf('Poles of G(s):\n');
disp(p);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
bode(num, den); grid on;
title('Bode Plot for Problem 1.1');

% Problem 1.2
fprintf('\nProblem 1.2\n');
num = [1.1057 -0.1900];
den = [1 17.95 123.3 366.3 112.2 0];
G = tf(num, den);
z = zero(G);
p = pole(G);
fprintf('Zeros of G(s):\n');
disp(z);
fprintf('Poles of G(s):\n');
disp(p);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
bode(num, den); grid on;
title('Bode Plot for Problem 1.2');

% Problem 2: B-7-3
den = [2 1];
% Part a
num = [1 1];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2: B-7-3, Part a');
% Part b
num = [-1 1];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 1: B-7-3, Part b');

% Problem 2: B-7-8
num = [-1 1]; den = [1 1];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2: B-7-8');

% Problem 2: B-7-10
num = [10 5];
den = [1 12 20 0 0];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2: B-7-10');

% Problem 2: B-7-12
num = 1;
den = [1 0.8 1 0];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2: B-7-12');

% Problem 2
T1 = 1/2; T2 = 1/10;    % set T1 > T2 > 0 to check (NOTE: T1/T2 = 5 here)
den = [T2 1];
% Part a
num = [T1 1];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2a');
% Part b
num = [T1 -1];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2b');
% Part c
num = [-T1 1];
G = tf(num, den);
figure('Units', 'inches', 'Position', [1 4 5 4]); 
nyquist(G); axis equal;
title('Nyquist Plot for Problem 2c');
