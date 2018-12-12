%%%%%
% Jordan Mayer
% AAE 564
% HW 14
%
% Perform various matrix operations to obtain and check results.
%%%%%

clear all; close all; format compact;

%%% Problem 5 %%%

P4 = [2 1 1 1 0.5 1];
E1 = [0 0 0];
E2 = [0 pi pi];

% set C, D such that measured output is y
C = [1 0 0 0 0 0];
D = 0;

fprintf('\nL7:\n\n');

[A_L7, B_L7] = get_A_B_double_pendulum(P4, E1)
if (is_controllable(A_L7, B_L7) && is_controllable(A_L7', C'))
    x_e_L7 = [E1, zeros(1,3)]'
    p_L7 = [-1 -1.2 -1.4 -1.6 -1.8 -2];
    K_L7 = -place(A_L7, B_L7, p_L7)
    L_L7 = -place(A_L7', C', p_L7)'
end

fprintf('\nL8:\n\n');

[A_L8, B_L8] = get_A_B_double_pendulum(P4, E2)
if (is_controllable(A_L8, B_L8) && is_controllable(A_L8', C'))
    x_e_L8 = [E2, zeros(1,3)]'
    p_L8 = [-1 -1.2 -1.4 -1.6 -1.8 -2]
    K_L8 = -place(A_L8, B_L8, p_L8)
    L_L8 = -place(A_L8', C', p_L8)'
end