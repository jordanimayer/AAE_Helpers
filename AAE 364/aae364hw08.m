%%%%%
% AAE 364
% HW 08
% Jordan Mayer
%
% Perform various calculations and create root loci for problems in
% HW 08.
% Output values (from disp) are written down in the hand-written
% portion of the assignment.
%%%%%

close all; % close all open figures

%%% Problem 1: B-6-12 %%%
% Part a
% Plot root locus in MATLAB to compare to hand-drawn sketch
L = tf([0 1 -1], [1 6 8]);  % L(s)
figure(1); rlocus(L);       % plot root locus
title('Root Locus for Problem 1: B-6-12, Part a');
axis([-6 6 -6 6]); pbaspect([1 1 1]);
    % square plot for easy comparison to hand-drawn sketch
    
% Part b
% Find break-in/break-away points of root locus
syms s;
L = (1-s)/((s+2)*(s+4));    % L(s)
D = -diff(1/L);
S = solve(D,s);             % solve -d/ds(1/L(s)) = 0
fprintf('\nProblem 1: B-6-12, Part b\nbreak-in/break-away points:\n');
disp(S);

% Plot root locus in MATLAB to compare to hand-drawn sketch
L = tf([0 -1 1], [1 6 8]);  % L(s)
figure(2); rlocus(L);       % plot root locus
title('Root Locus for Problem 1: B-6-12, Part b');
axis([-6 6 -6 6]); pbaspect([1 1 1]);
    % square plot for easy comparison to hand-drawn sketch

    
%%% Problem 2 %%%
% Part 1a
% Find zeroes and poles of L(s)
L = tf([0 0 1.1057 0.1900], [1 0.7385 0.8008 0]);   % L(s)
Z = zero(L); P = pole(L);
fprintf('\nProblem 2, Part 1a\nzeroes:\n');
disp(Z);
fprintf('Problem 2, Part 1a\npoles:\n');
disp(P);

% Part 1b
% Plot root locus in MATLAB to compare to hand-drawn sketch
figure(3); rlocus(L);   % plot root locus
title('Root Locus for Problem 2, Part 1b');
axis([-2 2 -2 2]); pbaspect([1 1 1]);
    % square plot for easy comparison to hand-drawn sketch

% Part 2a
% Find zeroes and poles of L(s)
L = tf([0 0 0 0 1.1057 -0.1900], [1 17.95 123.3 366.3 112.2 0]);    % L(s)
Z = zero(L); P = pole(L);
fprintf('\nProblem 2, Part 2a\nzeroes:\n');
disp(Z);
fprintf('Problem 2, Part 2a\npoles:\n');
disp(P);

% Find s such that -d/ds(1/L(s)) = 0
syms s;
N = poly2sym([0 0 0 0 1.1057 -0.1900], s);      % numerator of L(s)
D = poly2sym([1 17.95 123.3 366.3 112.2 0], s); % denominator of L(s)
Ls = N/D;   % L(s) as symbolic equation, not transfer function
D = diff(Ls);
S = vpasolve(-D == 0, s);
S = double(S);  % to avoid excessive decimal places
fprintf('Problem 2, Part 2a\nbreak-in/break-away points:\n');
disp(S);

% Solve system of equations for imaginary axis intersections
syms k w;
R = 17.95*w^4 - 366.3*w^2 - 0.1900*k;       % Real component of CE
I = w*(w^4 - 123.3*w^2 + (112.2+1.1057*k)); % Imaginary component of CE
[k, w] = solve([R == 0, I == 0], k, w);
k = double(k);  % to avoid excessive decimal places
w = double(w);  % to avoid excessive decimal places
fprintf('Problem 2, Part 2a\nimaginary axis intersection k:\n');
disp(k);
fprintf('\nProblem 2, Part 2a\nimaginary axis intersection w:\n');
disp(w);

% Part 2b
% Plot root locus in MATLAB to compare to hand-drawn sketch
figure(4); rlocus(L);   % plot root locus
title('Root Locus for Problem 2, Part 2b');
axis([-9 9 -9 9]); pbaspect([1 1 1]);
    % square plot for easy comparison to hand-drawn sketch