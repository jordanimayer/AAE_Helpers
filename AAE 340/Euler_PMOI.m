%%%%%
% AAE 340
% HW 10 - Problem 2
% Jordan Mayer
%
% Evaluate state space of differential equations for given system.
% Derived in hand-written portion.
%%%%%

function y_dot = Euler_PMOI(t, y)
    It = 2887;
    Ia = 5106;
    w_z = 1.1;
    lambda = (Ia-It)/It * w_z;

    % Obtain state variables
    y_1 = y(1);
    y_2 = y(2);
    y_3 = y(3);
    y_4 = y(4);
    y_5 = y(5);
    y_6 = y(6);
    
    % Calculate derivatives
    y_dot1 = y_4;
    y_dot2 = y_5;
    y_dot3 = y_6;
    y_dot4 = -lambda^2 * y_1;
    y_dot5 = -lambda^2 * y_2;
    y_dot6 = 0;
    
    % Construct derivative state vector
    y_dot = [y_dot1 y_dot2 y_dot3 y_dot4 y_dot5 y_dot6].';