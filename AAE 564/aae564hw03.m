%%%%%
% Jordan Mayer
% AAE 564
% HW 03
%
% Calculate state space realizations and check work
%%%%%

%%% Problem 2 and 3 %%%
% 2b: calculate A, B, C, D for each equilibrium state and parameters
% 3: calculate poles and zeroes of the system for each situation

% Initialize equilibrium states
u_e = 0;
y_e = 0;
theta1_e_list = [0 pi];
theta2_e_list = [0 pi];

% Initialize parameters [P1 P2 P3 P4]
m0_p = [2 2 2 2];
m1_p = [1 1 1 1];
m2_p = [1 1 0.5 1];
len1_p = [1 1 1 1];
len2_p = [1 0.99 1 0.5];
g_p = [1 1 1 1];

L = 0;  %  for text output only

% Iterate through parameters and equilibrium states
for P = [1 4]
    for E = [1 2]
        L = L + 1;
        fprintf('\nL%d:\n\n', L);
        
        % Get parameter and equilibrium values
        m0 = m0_p(P);
        m1 = m1_p(P);
        m2 = m2_p(P);
        len1 = len1_p(P);
        len2 = len2_p(P);
        g = g_p(P);
        theta1_e = theta1_e_list(E);
        theta2_e = theta2_e_list(E);
        
        a = m0 + (1 - (cos(theta1_e))^2) * m1 + ...
            (1 - (cos(theta2_e))^2) * m2;  % just a helper variable
        
        % calculate A, B, C, D for 2b
        A = [ ...
            0 1 0 0 0 0;
            0 0 -(1/a)*m1*g*(cos(theta1_e))^2 0 ...
              -(1/a)*m2*g*(cos(theta2_e))^2 0;
            0 0 0 1 0 0;
            0 0 -(1/(len1*a))*m1*g*(cos(theta1_e))^3 0 ...
              -(1/(len1*a))*m2*g*cos(theta1_e)*(cos(theta2_e))^2 0;
            0 0 0 0 0 1;
            0 0 -(1/(len2*a))*m1*g*(cos(theta1_e))^2*cos(theta2_e) 0 ...
              -(1/(len2*a))*m2*g*(cos(theta2_e))^3 0;
            ];
        B = [0 1/a 0 1/(len1*a)*cos(theta1_e) 0 1/(len2*a)*cos(theta2_e)]';
        C = [1 0 0 0 0 0];
        D = 0;
        
        % Display A, B, C, D for 2b
        fprintf('A:\n');
        disp(A);
        fprintf('B:\n');
        disp(B);
        fprintf('C:\n');
        disp(C);
        fprintf('D:\n');
        disp(D);
        % Note:
        % C and D will be the same for every case because the output is
        % always just y
        
        % Calculate and display poles, zeroes for 3
        [num, den] = ss2tf(A,B,C,D);
        G = tf(num,den);
        fprintf('poles:\n');
        disp(pole(G));
        fprintf('zeroes:\n');
        disp(zero(G));
        fprintf('_____\n');
    end
end