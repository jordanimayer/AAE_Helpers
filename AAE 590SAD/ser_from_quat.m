%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% ser_from_quat:
%   Calculate lambda, theta for single equivalent rotation (SRE) given a
%   quaternion between two frames.
%
% Inputs:
%   q: quaternion, [vector scalar]
%
% Outputs:
%   lambda: unit vector for axis of rotation, valid in either frame
%   theta: angle of rotation about lambda, deg
%%%%%

function [lambda, theta] = ser_from_quat(q)
    lambda = q(1:3)/norm(q(1:3));
    theta = 2*acosd(q(4));
end