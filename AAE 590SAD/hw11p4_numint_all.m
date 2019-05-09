%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 11
%
% hw11p4_sim_and_plot:
%   Calculate revolution rates of change for Euler parameters and
%   nondimensionalized angular velocities. Includes rotor (for use in
%   Problem 4).
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

function y_prime = hw11p4_numint_all(v, y)
    global I n_rotor;
    K = get_K_from_I(I);
    J = 82.5;  % rotor's axial moment of inerta, kg-met^2
    
    w1 = y(1);
    w2 = y(2);
    w3 = y(3);
    q1 = y(4);
    q2 = y(5);
    q3 = y(6);
    q4 = y(7);
    
    % calculate deREVatives (eh?)
    w1_prime = 2*pi*K(1)*(w2*w3 - ...
                          6*(1 - 2*q3^2 - 2*q1^2)*(q2*q3 - q1*q4));
    w2_prime = 2*pi*K(2)*(w3*w1 - ...
                          12*(q1*q2 + q3*q4)*(q2*q3 - q1*q4)) - ...
               2*pi*(J/I(2))*n_rotor*w3;
    w3_prime = 2*pi*K(3)*(w1*w2 - ...
                          6*(q1*q2 + q3*q4)*(1 - 2*q3^2 - 2*q1^2)) + ...
               2*pi*(J/I(3))*n_rotor*w2;
    q1_prime = pi*((w1 - 1)*q4 - w2*q3 + w3*q2);
    q2_prime = pi*((w1 + 1)*q3 + w2*q4 - w3*q1);
    q3_prime = pi*(-(w1 + 1)*q2 + w2*q1 + w3*q4);
    q4_prime = -pi*((w1 - 1)*q1 + w2*q2 + w3*q3);
    
    % fill state deREVative vector
    y_prime = [w1_prime w2_prime w3_prime ...
               q1_prime q2_prime q3_prime q4_prime].';
end