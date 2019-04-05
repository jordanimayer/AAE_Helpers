%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 05
%
% hw05_numint_qw:
%   Calculate time rates of change for both kinematic and dynamic
%   differential equations.
%
% Inputs:
%   t: current time, s
%   y: current state vector, [q1 q2 q3 q4 w1 w2 w3], where qi is ith
%      quaternion parameter and wi is ith angular velocity measure number,
%      with vectors expressed in s_hat frame
%
% Outputs:
%   y_dot: current state derivative vector (with respect to time)
%%%%%

function y_dot = hw05_numint_qw(t, y)
    % make typing a bit easier
    y1 = y(1);
    y2 = y(2);
    y3 = y(3);
    y4 = y(4);
    y5 = y(5);
    y6 = y(6);
    y7 = y(7);
    
    % set constants (described in handwritten work
    K1 = 0.7;
    K2 = 0.11;
    
    % actually calculate derivatives
    y_dot = zeros(7,1);
    y_dot(1) = 1/2 * (y5*y4 - y6*y3 + y7*y2);
    y_dot(2) = 1/2 * (y5*y3 + y6*y4 - y7*y1);
    y_dot(3) = 1/2 * (-y5*y2 + y6*y1 + y7*y4);
    y_dot(4) = -1/2 * (y5*y1 + y6*y2 + y7*y3);
    y_dot(5) = 0;
    y_dot(6) = K1*y7;
    y_dot(7) = -K1*y6 + K2;
end