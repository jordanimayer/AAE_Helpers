%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% rotvec_srt:
%   Express vector after simple rotation using Simple Rotation Theorem
%   (SRT). Result will be expressed in same reference frame as initial
%   (unrotated) vector.
%
% Inputs:
%   a: initial (unrotated) vector
%   lambda: unit vector along axis of rotation
%   theta: angle of rotation, deg
%
% Outputs:
%   b: rotated vector
%%%%%

function [b] = rotvec_srt(a, lambda, theta)
    b = a*cosd(theta) - cross(a,lambda)*sind(theta) + ...
        vec_dyadic_dot(a, dyadic_vecvec(lambda,lambda))*(1-cosd(theta));
end