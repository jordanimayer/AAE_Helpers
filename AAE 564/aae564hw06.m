%%%%%
% Jordan Mayer
% AAE 564
% HW 06
%
% Perform various linear algebra operations to confirm and support
% calculations
%%%%%

%%% Problem 1 %%%
fprintf('\nProblem 1:\n\n');
A = [2 -3 0; 2 -3 0; 3 -5 1];
% Use RREF to get eigenvectors
b = zeros(3,1);
for lambda = [0 1 -1]  % for each eigenvalue
    fprintf('lambda = %d\n', lambda);
    newA = lambda * eye(3) - A;
    fprintf('rref([newA b]) =\n');
    disp(rref([newA b]));
end

% Confirm eigenvalues and eigenvectors
[v, lambda] = eig(A)

%%% Problem 2 %%%
fprintf('\nProblem 2:\n\n');
A = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0];
% Use RREF to get eigenvectors
b = zeros(4,1);
for lambda = [1i -1i 1 -1]  % for each eigenvalue
    fprintf('lambda = %d + %dj\n', real(lambda), imag(lambda));
    newA = lambda * eye(4) - A;
    fprintf('rref([newA b]) =\n');
    disp(rref([newA b]));
end

% Confirm eigenvalues and eigenvectors
[v, lambda] = eig(A)

%%% Problem 3 %%%
fprintf('\nProblem 3:\n\n');
A = [0 1 0; 0 0 1; -1 1 1];
[v, lambda] = eig(A)
lambda1 = 1;  % eigenvalue with algebraic multiplicity 2
N = null(lambda1 * eye(3) - A);
nullity = size(N,2)

%%% Problem 4 %%%
p = [1 6 11 6];  % characteristic polynomial coefficients
A = compan(p)
[v, lambda] = eig(A)

%%% Problem 5 %%%
p = [1 -2 2];  % characteristic polynomial coefficients
A = compan(p)
[v, lambda] = eig(A)

%%% Problem 7 %%%
lambda = 2 + 3i;
mag = abs(lambda)
ang = angle(lambda)

%%% Problem 8 %%%

fprintf('\nProblem 8:\n\n');

% Parameters
m0 = [2 2 2 2];
m1 = [1 1 1 1];
m2 = [1 1 0.5 1];
len1 = [1 1 1 1];
len2 = [1 0.99 1 0.5];

% Equilibrium configurations
y_e = [0 0];
theta1_e = [0 pi];
theta2_e = [0 pi];

L = 0;
for p = 1:4
    for e = 1:2
        L = L + 1;
        fprintf('\nL%d:\n\n', L);
        m = [m0(p) m1(p) m2(p)]';
        len = [len1(p) len2(p)]';
        q_0 = [y_e(e) theta1_e(e) theta2_e(e)]';
        q_dot_0 = [0 0 0]';
        x_e = [q_0; q_dot_0];
        [A, B, C, D] = linmod('doublePendulum_nonlinear', x_e);
        fprintf('\nA =\n');
        disp(A);
        [v, lambda] = eig(A);
        for k = 1:6
            fprintf('lambda%d = %f + %fj\n', k,...
                    real(lambda(k,k)), imag(lambda(k,k)));
        end
    end
end

close all;  % close any open figures

m = [2 1 1]'; len = [1 0.5]';  % for nonlinear system
q_0 = rand(3,1); q_dot_0 = rand(3,1);
  % random initial state to compare results

% get nonlinear results
[t, x] = sim('doublePendulum_nonlinear');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 8b, nonlinear');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% get L7 results
[t, x] = sim('doublePendulum_L7');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 8b, L7');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% get L8 results
[t, x] = sim('doublePendulum_L8');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 8b, L8');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');


%%% Problem 9 %%%

% L7, stable solution
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 -0.5 -0.5 0 0 0;
     0 -1.5 -0.5 0 0 0;
     0 -1 -3 0 0 0];
  % matrix for L7 (see 8a)
[v, lambda] = eig(A);
x_0 = real(v(:,5));
q_0 = x_0(1:3);
q_dot_0 = x_0(4:6);

% get nonlinear results
[t, x] = sim('doublePendulum_nonlinear');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9a, nonlinear');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% get L7 results
[t, x] = sim('doublePendulum_L7');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9a, L7');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% L8, asymptotically stable solution
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 -0.5 -0.5 0 0 0;
     0 1.5 0.5 0 0 0;
     0 1 3 0 0 0];
  % matrix for L8 (see 8a)
[v, lambda] = eig(A);
x_0 = v(:,6);  % for some reason order is different than in
               % Problem 8 outputs
q_0 = x_0(1:3);
q_dot_0 = x_0(4:6);

% get nonlinear results
[t, x] = sim('doublePendulum_nonlinear');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9b, nonlinear');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% get L8 results
[t, x] = sim('doublePendulum_L8');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9b, L8');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% get L8 results for shorter amount of time
t = t(1:42);
y = y(1:42);
theta1 = theta1(1:42);
theta2 = theta2(1:42);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9b, L8 (shorter time)');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% L8, unstable solution
x_0 = v(:,4);  % for some reason order is different
q_0 = x_0(1:3);
q_dot_0 = x_0(4:6);

% get nonlinear results
[t, x] = sim('doublePendulum_nonlinear');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9c, nonlinear');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% get L8 results
[t, x] = sim('doublePendulum_L8');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9c, L8');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');

% get L8 results for shorter amount of time
t = t(1:42);
y = y(1:42);
theta1 = theta1(1:42);
theta2 = theta2(1:42);
figure();
subplot(3,1,1);
plot(t,y); grid on; xlabel('t'); ylabel('y');
title('Jordan Mayer - HW 06, Problem 9c, L8 (shorter time)');
subplot(3,1,2);
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
subplot(3,1,3);
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');
