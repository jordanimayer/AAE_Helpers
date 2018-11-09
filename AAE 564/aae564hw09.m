%%%%%
% Jordan Mayer
% AAE 564
% HW 09
%
% Perform basic matrix operations to obtain and check results
%%%%%

clear all; format compact;

%%% Problem 1 %%%

fprintf('\nProblem 1:\n');

fprintf('\nPart a:\n\n');

A = [-1 2001; -1 0];
[v, lambda] = eig(A)

fprintf('\nPart b:\n\n');

A = [-1 0; 0 1];
[v, lambda] = eig(A)

fprintf('\nPart c:\n\n');

A = [i 1; 0 i];
[v, lambda] = eig(A)
gm = size(null(lambda(1,1)*eye(2) - A),2)

fprintf('\nPart d:\n\n');
A = [-2 0; 0 0.5];
[v, lambda] = eig(A)

%%% Problem 2 %%%

fprintf('\nProblem 2:\n\n');
A = [1/2 1 -1/2 0; -1 1/2 0 -1/2; 1/2 0 -1/2 1; 0 1/2 -1 -1/2];
[v, lambda] = eig(A)
lambda1 = i;
gm1 = size(null(lambda1*eye(4) - A), 2)

%%% Problem 3 %%%

fprintf('\nProblem 3:\n');

fprintf('\nPart a:\n\n');
A = [0 1; -2 -sqrt(3)];
lambda = eig(A)

fprintf('\nPart b:\n\n');
A = [0 1; -1 0];
lambda = eig(A)

fprintf('\nPart c:\n\n');
A = [0 0 1 0; 0 0 0 1; 0 1 0 0; 1 0 0 0];
lambda = eig(A)

%%% Problem 4 %%%

fprintf('\nProblem 4:\n\n');
A = [0 1; -1/2 0];
lambda = eig(A)

%%% Problem 5 %%%

% set up parameters
m0 = 2; m1 = 1; m2 = 1; len1 = 1; len2 = 1; g = 1;

% set up equilibrium solutions
% NOTE: stability does not depend on y_e, as this does not affect A
theta1_e = [0 pi]; theta2_e = [0 pi];

% set up elements not dependent on equilibrium solution
M = zeros(3,3);
K = zeros(3,3);
M(1,1) = m0 + m1 + m2;
M(2,2) = m1 * len1^2;
M(3,3) = m2 * len2^2;

for k = [1 2]
    fprintf('L%d:\n\n', k);
    
    M(1,2) = -m1*len1*cos(theta1_e(k));
    M(1,3) = -m2*len2*cos(theta2_e(k));
    M(2,1) = M(1,2);
    M(3,1) = M(1,3);
    K(2,2) = -m1*len1*g*cos(theta1_e(k));
    K(3,3) = -m2*len2*g*cos(theta2_e(k));
    
    A = [zeros(3,3) eye(3); inv(M)*K zeros(3,3)]
    lambda = eig(A)
end

gm1 = size(null(A), 2)