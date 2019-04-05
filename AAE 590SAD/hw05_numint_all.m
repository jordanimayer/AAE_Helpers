%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 05
%
% hw05_numint_all:
%   Calculate time rates of change for Euler parameters, angular
%   velocities, AND direction cosine elements.
%
% Inputs:
%   t: current time, s
%   y: current state vector, [q1 q2 q3 q4 w1 w2 w3 C11 C21 C31]
%      where qi is the ith Euler parameter, with vector components
%      expressed in s_hat frame; wi is ith angular velocity measure number
%      (rad/s) as epxressed in s_hat frame; Cij is (i,j)th element of
%      direction cosine matrix from N to S.
%
% Outputs:
%   y_dot: current state derivative vector (with respect to time)
%%%%%

function y_dot = hw05_numint_all(t, y) 
    % extract state variables
    q1 = y(1);
    q2 = y(2);
    q3 = y(3);
    q4 = y(4);
    w1 = y(5);  % rad/s
    w2 = y(6);  % rad/s
    w3 = y(7);  % rad/s
    C = [y(8:10) y(11:13) y(14:16)].';
      % direction cosine matrix from N to S
      
    % set constants (described in handwritten work)
    K1 = 0.7;
    K2 = 0.11;
    
    % calculate derivatives
    q1_dot = 1/2 * (w1*q4 - w2*q3 + w3*q2);
    q2_dot = 1/2 * (w1*q3 + w2*q4 - w3*q1);
    q3_dot = 1/2 * (-w1*q2 + w2*q1 + w3*q4);
    q4_dot = -1/2 * (w1*q1 + w2*q2 + w3*q3);
    w1_dot = 0;
    w2_dot = K1*w3;
    w3_dot = -K1*w2 + K2;
    w_tilde = [0 -w3 w2; w3 0 -w1; -w2 w1 0];
    C_dot = C*w_tilde;
    
    % fill state derivative vector
    y_dot = [q1_dot q2_dot q3_dot q4_dot w1_dot w2_dot w3_dot ...
             C_dot(1,:) C_dot(2,:) C_dot(3,:)].';
end