%%%%%
% Jordan Mayer
% AAE 564
% HW 08
%
% Perform basic matrix arithmetic to obtain and check results
% Obtain linearizations and set up simulations
%%%%%

%%% Problem 1 %%%

fprintf('\nProblem 1:\n\n');

syms x;
C = [3 -1];
B = [0; 1];
xIminusA = [x -1; 2 x+3];
Ghat = C * inv(xIminusA) * B

%%% Problem 2 %%%

fprintf('\nProblem 2:\n\n');

C = [-1 -3];
B = [0; 1];
D = 1;
Ghat = C * inv(xIminusA) * B + D

%%% Problem 3 %%%

fprintf('\nProblem 3:\n\n');

% Part a

fprintf('\nPart a:\n\n');

A = [-5 2; -12 5];
[T, lambda] = eig(A)
B = [1 1]';
syms t;
eAt = expm(A*t)
x = expm(A*t) * B

% Part b

fprintf('\nPart b:\n\n');

B = [1 2]';
x = expm(A*t) * B

% Part c

fprintf('\nPart c:\n\n');

B = [1 3]';
x = expm(A*t) * B

%%% Problem 4 %%%

fprintf('\nProblem 4:\n\n');

A = [-1 0 0 0; 1 -2 0 0; 1 0 -3 0; 1 1 1 1];
B = [0 0 1 0]';
C = [0 0 1 0];

[T, lambda] = eig(A)
  % same values/vectors, but different order
lambda = diag([1 -1 -2 -3]);  % order in handwritten work
T = [T(:,1) T(:,4) T(:,2) T(:,3)];  % order in handwritten work
Thand = [0 4/sqrt(61) 0 0;...
         0 4/sqrt(61) 3/sqrt(10) 0;...
         0 2/sqrt(61) 0 4/sqrt(17);...
         1 -5/sqrt(61) -1/sqrt(10) -1/sqrt(17)];
  % handwritten T
Tdiff = T - Thand  % should be 0

invThand = [19/24 1/3 1/4 1;...
            sqrt(61)/4 0 0 0;...
            -sqrt(10)/3 sqrt(10)/3 0 0;...
            -sqrt(17)/8 0 sqrt(17)/4 0];
  % handwritten inv(T)
invTdiff = inv(T) - invThand  % should be 0

eAt = expm(A*t)
y = C * expm(A*t) * B

%%% Problem 5 %%%

fprintf('\nProblem 5:\n\n');

% Part a

fprintf('\nPart a:\n\n');

% Parameter set P4
m = [2 1 1];
len = [1 0.5];
g = 1;

% Equilibrium state E1
y_e = 0;
theta1_e = 0;
theta2_e = 0;
q_0 = [y_e theta1_e theta2_e]';
q_dot_0 = [0 0 0]';
x_e = [q_0; q_dot_0];

% Get linearization
k = 1;
c = 1;
alpha = 0;
omega = 1;
[A, B, C, D] = ...
  linmod('doublePendulum_nonlinear_passiveStabilizeDisturbance', x_e)
  % recall: u_e = 0
sys = ss(A,B,C,D);
p = pole(sys)
z = zero(sys)

% Part b

% choose disturbance such that steady-state response of linearized system
% is zero
eta = z(1);
alpha = real(eta);
omega = imag(eta);

% set up matrices for P4, E2, output y
B = [0 0 0 1/2 1/2 1]';
  % for P4, E1 (from HW 03 Solution)
C = [1 0 0 0 0 0];
D = 0;