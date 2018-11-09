%%%%%
% Jordan Mayer
% AAE 564
% HW 11
%
% Perform matrix operations to obtain and confirm results.
%%%%%

clear all; format compact;

%%% Problem 1 %%%

fprintf('\nProblem 1:\n');

% Part a

fprintf('\nPart a:\n\n');
A = [-1 0; 0 1]
C = [1 1]
Q0 = obsv(A,C)
r = rank(Q0)

% Part b

fprintf('\nPart b:\n\n');
C = [0 1]
Q0 = obsv(A,C)
r = rank(Q0)

% Part c

fprintf('\nPart c:\n\n');
A = [1 0; 0 1]
C = [1 1]
Q0 = obsv(A,C)
r = rank(Q0)

% Part d

fprintf('\nPart d:\n\n');
A = [0 1; 4 0]
C = [-2 1]
Q0 = obsv(A,C)
r = rank(Q0)

%%% Problem 2 %%%

fprintf('\nProblem 2:\n\n');

omega = 2; k = 1; m = 1;  % observability should not depend on these values
A = [0 1 0 0; omega^2 - k/m 0 k/m 0; 0 0 0 1; k/m 0 omega^2 - k/m 0]
C = [1 0 0 0]
Q0 = obsv(A,C)
r = rank(Q0)

%%% Problem 3 %%%

fprintf('\nProblem 3:\n');

% Part b

fprintf('\nPart b:\n\n');
Q0 = [0 1; 0 1]
N = null(Q0, 'r')

% Part c

fprintf('\nPart c:\n\n');
Q0 = [1 1; 1 1]
N = null(Q0, 'r')

%%% Problem 4 %%%

fprintf('\nProblem 4:\n')

% Part b

fprintf('\nPart b:\n\n')
A = [-1 0; 0 1]
C = [0 1]
lambda = eig(A)
M1 = [A - lambda(1)*eye(2); C]
r1 = rank(M1)
M2 = [A - lambda(2)*eye(2); C]
r2 = rank(M2)

% Part c

fprintf('\nPart c:\n\n');

A = eye(2)
C = [1 1]
lambda = eig(A)
M = [A - lambda(1)*eye(2); C]
r = rank(M)

% Part d

fprintf('\nPart d:\n\n')

A = [0 1; 4 0]
C = [-2 1]
lambda = eig(A)
M1 = [A - lambda(1)*eye(2); C]
r1 = rank(M1)
M2 = [A - lambda(2)*eye(2); C]
r2 = rank(M2)

%%% Problem 5 %%%

fprintf('\nProblem 5:\n\n');

A = [5 -1 -2; 1 3 -2; -1 -1 4]
C = [1 1 0; 0 1 1]
Q0 = obsv(A,C)
r = rank(Q0)

lambda = eig(A)
lambda1 = lambda(3)
lambda2 = lambda(2)
lambda3 = lambda(1)
M1 = [A - lambda1*eye(3); C]
r1 = rank(M1)
M2 = [A - lambda2*eye(3); C]
r2 = rank(M2)
M3 = [A - lambda3*eye(3); C]
r3 = rank(M3)

%%% Problem 7 %%%

fprintf('\nProblem 7:\n');

P1 = [2 1 1 1 1 1];
P2 = [2 1 1 1 0.99 1];
P4 = [2 1 1 1 0.5 1];
E1 = [0 0 0];

C = [1 0 0 0 0 0];

% L1

fprintf('\nL1:\n');

A = get_A_double_pendulum(P1, E1)
Q0 = obsv(A,C)
rq = rank(Q0)
sigma = svd(Q0)
lambda = eig(A)

for k = 1:6
    fprintf('lambda(%d):\n', k);
    M = [A - lambda(k)*eye(6); C];
    r = rank(M)
end

% L3

fprintf('\nL3:\n');

A = get_A_double_pendulum(P2, E1)
Q0 = obsv(A,C)
rq = rank(Q0)
sigma = svd(Q0)

% L7

fprintf('\nL3:\n');

A = get_A_double_pendulum(P4, E1)
Q0 = obsv(A,C)
rq = rank(Q0)
sigma = svd(Q0)

%%% Problem 8 %%%

fprintf('\nProblem 8:\n\n');

% set k, m, omega such that sqrt(k/(2*m)) > omega
k = 4
m = 0.5
omega = sqrt(k/(2*m)) - 1.5
A = [0 1 0 0; omega^2 - k/m 0 k/m 0; 0 0 0 1; k/m 0 omega^2 - k/m 0]
C = [1/2 0 1/2 0]

% Part a

fprintf('\nPart a:\n\n');
Q0 = obsv(A,C)
N = null(Q0, 'r')

lambda = eig(A)
lambda1 = omega;
lambda2 = -omega;
lambda3 = lambda(1);
lambda4 = lambda(2);
r1 = rank([A - lambda1*eye(4); C])
r2 = rank([A - lambda2*eye(4); C])
r3 = rank([A - lambda3*eye(4); C])
r4 = rank([A - lambda4*eye(4); C])