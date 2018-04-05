%%%%%
% AAE 340 - HW 07
% Jordan Mayer
%%%%%

%%% Problem 2 %%%
% Evaluate state space system of differential equations for binary star
% system. Described in 2b.

function y_dot = binarystar (t, y)
    % Initialize known values
    G = 6.673*10^(-11);         % universal gravitational constant, m^3/(kg*s)
    m_sun = 1.987*10^30;        % mass of sun, kg
    % Calculate m1, m2
    m1 = m_sun;                 % mass of star 1, kg
    m2 = m1/2;                  % mass of star 2, kg
    
    % Obtain state variables
    y_1 = y(1);
    y_2 = y(2);
    y_3 = y(3);
    y_4 = y(4);
    y_5 = y(5);
    y_6 = y(6);
    y_7 = y(7);
    y_8 = y(8);
    
    % Calculate derivative of state variables
    y_dot1 = y_5;
    y_dot2 = y_6;
    y_dot3 = y_7;
    y_dot4 = y_8;
    y_dot5 = G*m2/((y_3-y_1)^2 + (y_4-y_2)^2)^(3/2)*(y_3-y_1);
    y_dot6 = G*m2/((y_3-y_1)^2 + (y_4-y_2)^2)^(3/2)*(y_4-y_2);
    y_dot7 = G*m1/((y_3-y_1)^2 + (y_4-y_2)^2)^(3/2)*(y_1-y_3);
    y_dot8 = G*m1/((y_3-y_1)^2 + (y_4-y_2)^2)^(3/2)*(y_2-y_4);
    
    % Put state variables into derivative state vector
    y_dot = [y_dot1 y_dot2 y_dot3 y_dot4 y_dot5 y_dot6 y_dot7 y_dot8]';
end