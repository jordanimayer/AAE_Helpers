%%%%%
% Jordan Mayer
% AAE 364
% HW 10
%
% Create bode plots to confirm hand-drawn sketches.
%%%%%

close all;  % close any open figures

% Problem 1: B-7-3
% Part a
num = [1 1]; den = [2 1];
figure(); bode(num, den);
title('Bode Plot for Problem 1: B-7-3, Part a'); grid on;
% Part b
num = [-1 1]; den = [2 1];
figure(); bode(num, den);
title('Bode Plot for Problem 1: B-7-3, Part b'); grid on;

% Problem 1: B-7-8
num = [-1 1]; den = [1 1];
figure(); bode(num, den);
title('Bode Plot for Problem 1: B-7-8'); grid on;

% Problem 1: B-7-10
num = [10 5]; den = [1 12 20 0 0];
figure(); bode(num, den);
title('Bode Plot for Problem 1: B-7-10'); grid on;

% Problem 1: B-7-12
num = [0 0 0 1]; den = [1 0.8 1 0];
figure(); bode(num, den);
title('Bode Plot for Problem 1: B-7-12'); grid on;

% Problem 2.1
T1 = 1/10; T2 = 1/100;    % set T1 > T2 > 0 to check
% Part a
num = [T1 1]; den = [T2 1];
figure(); bode(num, den);
title('Bode Plot for Problem 2.1, Part a'); grid on;
% Part b
num = [T1 -1]; den = [T2 1];
figure(); bode(num, den);
title('Bode Plot for Problem 2.1, Part b'); grid on;
% Part c
num = [-T1 1]; den = [T2 1];
figure(); bode(num, den);
title('Bode Plot for Problem 2.1, Part c'); grid on;

% Problem 2.2
% Part a
Ta = 1/10; Tb = 1/100; T = 1/1000;  % set Ta > Tb > T > 0 to check
    % set K = 1 to check
num = [Ta*Tb Ta+Tb 1]; den = [T 1 0 0];
figure(); bode(num, den);
title('Bode Plot for Problem 2.2, Part a'); grid on;
% Part b
T = 1/10; Ta = 1/100; Tb = 1/1000;  % set T > Ta > Tb > 0 to check
num = [Ta*Tb Ta+Tb 1]; den = [T 1 0 0];
figure(); bode(num, den);
title('Bode Plot for Problem 2.2, Part b'); grid on;
% Part c
Ta = 1/10; T = 1/100; Tb = 1/1000;  % set Ta > T > Tb > 0 to check
num = [Ta*Tb Ta+Tb 1]; den = [T 1 0 0];
figure(); bode(num, den);
title('Bode Plot for Problem 2.2, Part c'); grid on;