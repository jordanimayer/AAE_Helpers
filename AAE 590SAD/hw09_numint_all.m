%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 09
%
% hw09_numint_all:
%   Calculate revolution rates of range for Euler parameters and
%   nondimensionalized angular velocities.
%
% Inputs:
%   v: current number of revolutions
%   y: current state vector, [w1 w2 w3 q1 q2 q3 q4]
%      wi is the ith nondimensionzlied angular velocity of S in N, with
%      vector components expressed in the c_hat frame
%      qi is the ith Euler parameter from the A to the C frame, with
%      vector components expressed in the c_hat frame
%
% Outputs:
%   y_prime: current state derivative vector (with respect to number of
%            revolutions)
%%%%%

function y_prime = hw09_numint_all(v, y)
    % extract state variables
    w1 = y(1);
    w2 = y(2);
    w3 = y(3);
    q1 = y(4);
    q2 = y(5);
    q3 = y(6);
    q4 = y(7);
    
    % calculate extra variables
    J = 150;  % principal moment of inertia about body axis of symmetry,
              % kg*m^2
    I = 500;  % principal moment of inertia about transverse body axes,
              % kg*m^2
    x = (J-I)/I;
    y = w1 - 1;  % note that w1 should never change since w1_prime = 0
    
    % calculate deREVatives (eh?)
    w1_prime = 0;
    w2_prime = -2*pi*y*w3 + ...
               2*pi*x*(12*(q1*q2 + q3*q4)*(q2*q3 - q1*q4) - w1*w3);
    w3_prime = 2*pi*y*w2 - ...
               2*pi*x*(6*(q1*q2 + q3*q4)*(1 - 2*q3^2 - 2*q1^2) - w1*w2);
    q1_prime = pi*(q4*(w1 - y - 1) - q3*w2 + q2*w3);
    q2_prime = pi*(q3*(w1 - y + 1) + q4*w2 - q1*w3);
    q3_prime = pi*(-q2*(w1 - y + 1) + q1*w2 + q4*w3);
    q4_prime = -pi*(q1*(w1 - y -1) + q2*w2 + q3*w3);
    
    % fill state deREVative vector
    y_prime = [w1_prime w2_prime w3_prime ...
               q1_prime q2_prime q3_prime q4_prime].';
end