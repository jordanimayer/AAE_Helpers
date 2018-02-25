% AAE 301 - HW 09
%
% Jordan Mayer

% Page 220-222
% Problem 1

% Part 1
M = [1 0 0; 0 2 0; 0 0 3];
phi = [6 0 -1; 0 0 0; -1 0 1];
K = [4 -1 0; -1 5 -2; 0 -2 2];
b = [0; 1; 0];

% Part 2
disp(M == conj(M))
eigM = eig(M)
disp(K == conj(K))
eigK = eig(K)
disp(phi == conj(phi))
eigphi = eig(phi)

% Part 3
A = [zeros(3,3), eye(3); -inv(M)*K, -inv(M)*phi]
B = [zeros(3,1); inv(M)*b]

% Part 4
syms s;
char = det(s*eye(6) - A)
eig(A)

% Part 5
u = 3;
q_inf = inv(K)*b*u
[num, den] = ss2tf(A, B, [eye(3), zeros(3,3)], zeros(3,1));

% Parts 5 and 6
G = [minreal(tf(num(1,:),den)); minreal(tf(num(2,:),den)); minreal(tf(num(3,:),den))]
step(u*G(1));hold on;step(u*G(2));step(u*G(3));
title('Step response for qdot when u(t) = 3 (pg. 220 #1)');
legend('q1', 'q2', 'q3');

% Part 7
[num, den] = ss2tf(A,B,-[K,phi],b);
H = [minreal(tf(num(1,:),den)); minreal(tf(num(2,:),den)); minreal(tf(num(3,:),den))]
