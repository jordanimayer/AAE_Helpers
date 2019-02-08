%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% srd_lambdatheta:
%   Calculate simple rotation dyadic from dextral orthonormal vector basis
%   a to basis b due to simple rotation about axis lambda by angle theta.
%   Assumes b-frame was equal to a-frame prior to rotation.
%
% Inputs:
%   lambda: rotation axis unit vector, in a-frame (and b-frame)
%           [lambda1 lambda2 lambda3]
%   theta: rotation angle, deg
%
% Outputs:
%   R: simple rotation dyadic (SRD) from a-frame to b-frame, matrix form
%%%%%

function [R] = srd_lambdatheta(lambda, theta)
    R = eye(3)*cosd(theta) - ...
        dyadic_vec_cross(eye(3), lambda)*sind(theta) + ...
        dyadic_vecvec(lambda,lambda)*(1 - cosd(theta));
end