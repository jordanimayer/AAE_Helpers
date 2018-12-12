%%%%%
% Jordan Mayer
% AAE 564
% HW 12
%
% Perform matrix operations to obtain and confirm results.
%%%%%

clear all; format compact; rehash;

fprintf('\nProblem 1:\n');

fprintf('\nPart a\n\n');

A = [-1 0; 0 1]
B = [1 1]'
Qc = ctrb(A, B)
r = rank(Qc)

fprintf('\nPart b:\n\n');

A = [-1 0; 0 1]
B = [0 1]'
Qc = ctrb(A, B)
r = rank(Qc)

fprintf('\nPart c:\n\n');

A = [1 0; 0 1]
B = [1 1]'
Qc = ctrb(A, B)
r = rank(Qc)

fprintf('\nProblem 2:\n\n');

A = [5 1 -1; -1 3 -1; -2 -2 4]
B = [1 0; 1 1; 0 1]
Qc = ctrb(A, B)
r = rank(Qc)
lambda = eig(A)
r2 = rank([A - lambda(1)*eye(3), B])
r6 = rank([A - lambda(2)*eye(3), B])
r4 = rank([A - lambda(3)*eye(3), B])

fprintf('\nProblem 3:\n');

P1 = [2 1 1 1 1 1];  % from HW 01
P2 = [2 1 1 1 0.99 1];  % from HW 01
P4 = [2 1 1 1 0.5 1];  % from HW 01
E1 = [0 0 0];  % from HW 06
phi = [1 0 0];  % from HW 06 solution

fprintf('\nL1:\n\n');

[A, B] = get_A_B_double_pendulum(P1, E1)
Qc = ctrb(A, B)
r = rank(Qc)
sigma = svd(Qc)
lambdas = eig(A)
for k = 1:6
    fprintf('\n');
    lambda = lambdas(k)
    r_pbh = rank([A - lambda*eye(6), B])
end

fprintf('\nL3:\n\n')

[A, B] = get_A_B_double_pendulum(P2, E1)
Qc = ctrb(A, B)
r = rank(Qc)
sigma = svd(Qc)
lambdas = eig(A)
for k = 1:6
    fprintf('\n');
    lambda = lambdas(k)
    r_pbh = rank([A - lambda*eye(6), B])
end

fprintf('\nL7:\n\n')

[A, B] = get_A_B_double_pendulum(P4, E1)
Qc = ctrb(A, B)
r = rank(Qc)
sigma = svd(Qc)
lambdas = eig(A)
for k = 1:6
    fprintf('\n');
    lambda = lambdas(k)
    r_pbh = rank([A - lambda*eye(6), B])
end

fprintf('\nProblem 4:\n\n');

omega = rand; k = rand; m = rand;
A = [0 1 0 0; omega^2 - k/m 0 k/m 0; 0 0 0 1; k/m 0 omega^2 - k/m 0]
B = [0 0 0 1]'
Qc = ctrb(A, B)
r = rank(Qc)

fprintf('\nProblem 5:\n\n');

k = rand 
m = rand
omega = sqrt(k/(2*m)) - rand
A = [0 1 0 0; omega^2 - k/m 0 k/m 0; 0 0 0 1; k/m 0 omega^2 - k/m 0]
B = [0 -1 0 1]'
Qc = ctrb(A, B)
r = rank(Qc)
lambda = eig(A)
lambda1 = omega
lambda2 = -omega
lambda3 = sqrt(2*k/m - omega^2) * 1i
lambda4 = -sqrt(2*k/m - omega^2) * 1i
r1 = rank([A - lambda1*eye(4), B])
r2 = rank([A - lambda2*eye(4), B])
r3 = rank([A - lambda3*eye(4), B])
r4 = rank([A - lambda4*eye(4), B])

fprintf('\n');

% check basis for controllable subspace
% (should have same range as orth(Qc))
t1t2 = [0 -1 0 1; -1 0 1 0]'
r1 = rank(t1t2)
r2 = rank(orth(Qc))
r = rank([t1t2 orth(Qc)])  % same as r1 and r2 (we're good!)

T = [0 -1 0 1; -1 0 1 0; 0 1 0 1; 1 0 1 0]
rT = rank(T)  % should be 4 (it is!)
Tinv = inv(T)

syms w k m;  % w represents omega
A = [0 1 0 0; w^2 - k/m 0 k/m 0; 0 0 0 1; k/m 0 w^2 - k/m 0]
B = [0 -1 0 1]'
C = [1 0 0 0]
A_tilde = inv(T) * A * T
B_tilde = inv(T) * B
C_tilde = C * T

% test controllable subsystem w/ random omega, k, m
k = rand
m = rand
w = sqrt(k/(2*m)) - rand
A = [0 1 0 0; w^2 - k/m 0 k/m 0; 0 0 0 1; k/m 0 w^2 - k/m 0]
B = [0 -1 0 1]'
C = [1 0 0 0]
Acc = [0 w^2 - 2*k/m; 1 0]
Bc = [1 0]'
Cc = [0 -1]
% controllable?
rc = rank(ctrb(Acc, Bc))  % should be 2 (it is!)
% same I-O behavior?
sys = ss(A, B, C, 0);
z = zero(sys)
p = pole(sys)  % 2 of these should be zeros also (they are!)
sysc = ss(Acc, Bc, Cc, 0);
zc = zero(sysc)
pc = pole(sysc)  % these should be the non-zero poles from sys (they are!)
% controllable eigenvalues?
hw_lambdac = sqrt(2 * k/m - w^2) * 1i  % from handwritten work
lambdac = eig(Acc)  % should be plus/minus hw_lambdac (they are!)

fprintf('\nProblem 7:\n\n');

A = [0 1; 4 0]
[v, d, w] = eig(A)