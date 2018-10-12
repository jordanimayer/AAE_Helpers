%%%%%
% Jordan Mayer
% AAE 564
% HW 05
% Problem 2
%
% Check null space and nullity of given matrix, A
%%%%%

A = [1 2 3 4; 2 3 4 5; 3 4 5 6];
basis_nullSpace_A = null(A, 'r')  % 'r' requests rational basis, rather
                                  % than orthonormal, which is better for
                                  % checking our handwritten work
nullity_A = size(basis_nullSpace_A, 2)