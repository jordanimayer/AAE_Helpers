% AAE 301 - HW 07
% Jordan Mayer
% Page 147

% Problem 2
% A = [0, -5; 1, -2];
% Problem 3
% A = [0, 1, 0; 0, 0, 1; 0, -4, 0]
% Problem 5
% A = [5, 8; -5, -7]

% Problem 6
A = [0, -sqrt(2), 0; sqrt(2), 0, -sqrt(2); 0, sqrt(2), 0];
syms t; syms s;
si = s*eye(3,3) - A
dsi = det(si)
isi = inv(si)
simplify(expm(A*t))