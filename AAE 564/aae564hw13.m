%%%%%
% Jordan Mayer
% AAE 564
% HW 13
%
% Perform various math operations to obtain and verify results.
%%%%%

clear all; close all; format compact;

%%% Problem 5 %%%

P2 = [2 1 1 1 0.99 1];
P4 = [2 1 1 1 0.5 1];
E1 = [0 0 0];
E2 = [0 pi pi];

% set C, D so y = x for simulation
C = eye(6);
D = zeros(6,1);

fprintf('\nL3:\n\n');

[A_L3, B_L3] = get_A_B_double_pendulum(P2, E1)
if (is_controllable(A_L3, B_L3))
    p_L3 = [-0.1067   -0.9619   -0.0046   -0.7749   -0.8173   -0.8687];
    K_L3 = -place(A_L3, B_L3, p_L3)
    x_e_L3 = [E1, zeros(1,3)]'
end

fprintf('\nL4:\n\n');

[A_L4, B_L4] = get_A_B_double_pendulum(P2, E2)
if (is_controllable(A_L4, B_L4))
    p_L4 = [-0.2785   -0.5469   -0.9575   -0.9649   -0.1576   -0.9706];
    K_L4 = -place(A_L4, B_L4, p_L4)
    x_e_L4 = [E2, zeros(1,3)]'
end

fprintf('\nL7:\n\n');

[A_L7, B_L7] = get_A_B_double_pendulum(P4, E1)
if (is_controllable(A_L7, B_L7))
    p_L7 = [-0.6948   -0.3171   -0.9502   -0.0344   -0.4387   -0.3816];
    K_L7 = -place(A_L7, B_L7, p_L7)
    x_e_L7 = [E1, zeros(1,3)]'
end

fprintf('\nL8:\n\n');

[A_L8, B_L8] = get_A_B_double_pendulum(P4, E2)
if (is_controllable(A_L8, B_L8))
    p_L8 = [-0.7655   -0.7952   -0.1869   -0.4898   -0.4456   -0.6463]
    K_L8 = -place(A_L8, B_L8, p_L8)
    x_e_L8 = [E2, zeros(1,3)]'
end

%%% Problem 7 %%%

fprintf('\nProblem 7:\n\n');

A = [0 1 0 0; -1/3 0 2/3 0; 0 0 0 1; 2/3 0 -1/3 0];
B = [0 0 ; 2/3 -1/3; 0 0; -1/3 2/3];
r = rank(ctrb(A,B))
p = rand(1,4) * -1  % eigenvalues all negative real
K = -place(A,B,p)

%%% Problem 8 %%%
Ac = [1 0; 0 -1];
Bc = [1 2]';
Kd = [-8/3 1/6];
T = log(2);

%%% Problem 10 %%%

A = [0 -1; -1 0];
B = [1 -1]';
C = [1 -1];
D = 0;
L = [1 3]';