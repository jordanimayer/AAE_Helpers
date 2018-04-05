%%%%%
% AAE 340
% HW 08
% Jordan Mayer
%%%%%

%%% Problem 1 %%%
% Evaluate state space system of differential equations for given system.
% Described in part a (hand-written).

function y_dot = block_rod_system(t, y)
    % Initialize known values
    g = 9.81;   % graviational constant, m/s^2
    L = 5;      % length of rod, m
    
    % Obtain state variables
    y_1 = y(1);
    y_2 = y(2);
    y_3 = y(3);
    y_4 = y(4);
    
    % Calculate derivative of state variables
    y_dot1 = y_3;
    y_dot2 = y_4;
    y_dot4 = (-g/L*sin(y_2)-sin(y_2)*cos(y_2)/2*y_4^2)...
        /(1-(cos(y_2))^2/2);
    y_dot3 = L/2*y_4^2*sin(y_2) - L/2*y_dot4*cos(y_2);
    
    % Put state variables into derivative state vector
    y_dot = [y_dot1 y_dot2 y_dot3 y_dot4].';
end