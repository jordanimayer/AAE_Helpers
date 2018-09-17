%%%%%
% AAE 564
% HW 01
% Jordan Mayer
%%%%%

%%% Problem 7 %%%
% Evaluate state space system of differential equations for
% double-pendulum system.

function x_dot = double_pendulum(t, x)
    P = 4;  % CHANGE THIS to change parameter set

    % Initialize parameters [P1 P2 P3 P4]
    m0_p = [2 2 2 2];
    m1_p = [1 1 1 1];
    m2_p = [1 1 0.5 1];
    len1_p = [1 1 1 1];
    len2_p = [1 0.99 1 0.5];
    g_p = [1 1 1 1];
    m0 = m0_p(P);
    m1 = m1_p(P);
    m2 = m2_p(P);
    len1 = len1_p(P);
    len2 = len2_p(P);
    g = g_p(P);
    
    % Set input
    u = 0;
    % Obtain state variables
    x_1 = x(1);
    x_2 = x(2);
    x_3 = x(3);
    x_4 = x(4);
    x_5 = x(5);
    x_6 = x(6);
    
    % Calculate state derivatives
    x_dot1 = x_2;
    x_dot2 = (-m1*sind(x_3) * (g*cosd(x_3) + len1*x_4^2) -...
              m2*sind(x_5) * (g*cosd(x_5) + len2*x_6^2) + u) /...
              (m0 + m1*(1 - (cosd(x_3))^2) + m2*(1 - (cosd(x_5))^2));
    x_dot3 = x_4;
    x_dot4 = (1/len1) * (cosd(x_3)*x_dot2 - g*sind(x_3));
    x_dot5 = x_6;
    x_dot6 = (1/len2) * (cosd(x_5)*x_dot2 - g*sind(x_5));
    
    % Put state derivatives into state derivative vector
    x_dot = [x_dot1 x_dot2 x_dot3 x_dot4 x_dot5 x_dot6].';
end