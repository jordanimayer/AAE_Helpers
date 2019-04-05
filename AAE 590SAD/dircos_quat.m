%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% dircos_quat:
%   Obtain direction cosine matrix using Euler parameters. Note that Euler
%   parameter vector component must be expressed in terms of either frame A
%   or frame B, where the direction cosine matrix is from frame A to frame
%   B.
%
% Inputs:
%   q: Euler parameter, [vector, scalar]
%
% Outputs:
%   C: direction cosine matrix
%%%%%

function [C] = dircos_quat(q)
    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    C = zeros(3,3);
    
    C(1,1) = 1 - 2*q2^2 - 2*q3^2;
    C(1,2) = 2*(q1*q2 - q3*q4);
    C(1,3) = 2*(q3*q1 + q2*q4);
    C(2,1) = 2*(q1*q2 + q3*q4);
    C(2,2) = 1 - 2*q3^2 - 2*q1^2;
    C(2,3) = 2*(q2*q3 - q1*q4);
    C(3,1) = 2*(q3*q1 - q2*q4);
    C(3,2) = 2*(q2*q3 + q1*q4);
    C(3,3) = 1 - 2*q1^2 - 2*q2^2;
end