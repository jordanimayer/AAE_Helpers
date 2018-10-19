%%%%%
% Jordan Mayer
% AAE 564
% HW 07
%
% Perform basic matrix arithmetic to obtain and check results
%%%%%

%%% Problem 2 %%%

fprintf('\nProblem 2:\n\n');

A = [0 1 0; 0 0 1; -1 -1 -2];  % companion matrix fitting description
I = eye(3);
fprintf('Check A^3\n');
disp(A^3 - (-I - A - 2*A^2));  % should be 0
fprintf('Check A^5\n');
disp(A^5 - (-3*I - A - 5*A^2));  % should be 0
fprintf('Check A^(-1)\n');
disp(inv(A) - (-A^2 - 2*A - I));  % should be 0

%%% Problem 3 %%%

fprintf('\nProblem 3:\n\n');

A = [0 1 0 0; 0 0 1 0; 0 0 0 1; 564 1 0 0];
fprintf('Check A^4\n');
disp(A^4);

%%% Problem 5 %%%

fprintf('\nProblem 5:\n\n');

A = [0 1; 1 0];
fprintf('e^(Aln2)\n');
disp(expm(A*log(2)));

%%% Problem 6 %%%

fprintf('\nProblem 6:\n\n');

A = [3 -5; -5 3];
fprintf('e^(At)\n');
syms t;
disp([-1 1; 1 1] * [exp(8*t) 0; 0 exp(-2*t)] * inv([-1 1; 1 1]));