%%%%%
% Jordan Mayer
% AAE 564
% HW 10
%
% Perform matrix arithmetic to check results.
%%%%%

clear all; format compact;

%%% Problem 1 %%%

fprintf('\nProblem 1:\n\n');

A = [1 j 0; -j 2 1; 0 1 4]
lambda = eig(A)

A = [1 j; -j 1]
lambda = eig(A)

A = [0 2; 2 0]
lambda = eig(A)

A = [-1 1; 1 -2]
negLambda = eig(-A)

fprintf('\nProblem 2:\n\n');

A = [3 4]'
sigma = svd(A)

A = [3 4]
sigma = svd(A)

A = [0 1; 1 0]
sigma = svd(A)

%%% Problem 3 %%%

fprintf('\nProblem 3:\n\n');

A = [3 0 1; 1 0 3]
U = 1/sqrt(2) * [1 1; -1 1]
S = [2 0 0; 0 4 0]
V = 1/sqrt(2) * [1 1 0; 0 0 sqrt(2); -1 1 0]
svdA = U * S * V'