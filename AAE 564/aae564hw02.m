%%%%%
% Jordan Mayer
% AAE 564
% HW 02
%
% Perform matrix arithmetic to compute transfer functions.
% See handwritten work for previous steps.
%%%%%


%%% Problem 5 %%%
A = [0 1 0 0; -2 0 -4 0; 0 0 0 1; 1 0 3 0];
B = [0 0 0 0; 3 0 0 3; 0 0 0 0; -2 0 0 -2];
C = [1 0 0 0; 0 0 1 0];
D = zeros(2,4);

syms s;
G = C * inv(s * eye(4) - A) * B + D;
fprintf('\nProblem 5:\nG =\n');
disp(G);

%%% Problem 6 %%%
A = [0 1 0; 0 -4 1; 0 1 -1];
B = [ 0 0; -5 4; -3 -1];
C = [1 0 1];
D = [0 0];

G = C * inv(s * eye(3) - A) * B + D;
fprintf('\nProblem 6:\nG =\n');
disp(G);