%%%%%
% Jordan Mayer
% AAE 564
% HW 05
% Problem 3
%
% Check range and rank of given matrix, A
%%%%%

A = [1 2 3 4; 2 3 4 5; 3 4 5 6];
basis_range_A = orth(A)
rank_A = rank(A)

B = [1 2; 2 3; 3 4]  % hand-calculated basis for range(A)
B_tilde = basis_range_A  % MATLAB-calculated orthonormal basis for range(A)
rank_B = rank(B)
rank_B_tilde = rank(B_tilde)
rank_B_B_tilde = rank([B B_tilde])