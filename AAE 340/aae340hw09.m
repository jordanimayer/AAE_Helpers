%%%%%
% AAE 340
% HW 09
% Jordan Mayer
%
% Perform matrix arithmetic to find inertia matrices.
%%%%%

%%% Problem 1a %%%
% compute inertia matrix in s-frame about point A
I_u_A = [216 -54 0; -54 18 0; 0 0 234];
    % inertia matrix in u-frame about point A
theta = atan(12/6); % angle between u1 and s1, rad
l_u_s = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
    % rotation matrix from u-frame to s-frame
I_s_A = l_u_s*I_u_A*l_u_s.';
    % inertia matrix in s-frame about point A
fprintf('\nProblem 1a\ninertia matrix in s-frame about point A\n');
disp(I_s_A);

%%% Problem 1d %%%
% factor equation and check eigenvalues/eigenvectors
syms x;     % variable for eigenvalues
d = (14.40-x)*(111.7-x)*(126.1-x)-0+0-25.22^2*(126.1-x)+0-0;
    % expression to solve for eigenvalues
f = factor(d, 'FactorMode', 'real');  % factors of expression
fprintf('\nProblem 1d\nfactors of eigenvalue expression\n');
disp(f);
I_s_B = [14.40 25.22 0; 25.22 111.7 0; 0 0 126.1];
    % inertia matrix in s-frame relative to point B
[vecs, vals] = eig(I_s_B);  % [eigenvalues, eigenvectors]
fprintf('\neigenvalues of I_s_B\n');
disp(vals);
fprintf('\neigenvectors of I_s_B\n');
disp(vecs);

%%% Problem 1f %%%
I_dash = vecs*I_s_B*vecs.';
fprintf('\nProblem 1f\ninertia matrix in e-frame about point B\n');
disp(I_dash);