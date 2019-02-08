%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% dircos_lambdatheta:
%   Calculate direction cosine matrix from dextral orthonormal vector basis
%   a to basis b due to simple rotation about axis lambda by angle theta.
%   Assumes b-frame was equal to a-frame prior to rotation.
%
% Inputs:
%   lambda: rotation axis unit vector, in a-frame (and b-frame)
%           [lambda1 lambda2 lambda3]
%   theta: rotation angle, deg
%
% Outputs:
%   aCb: direction cosine matrix from a-frame to b-frame
%%%%%

function [aCb] = dircos_lambdatheta(lambda, theta)
    aCb = zeros(3,3);
    
    for r = 1:3  % row of aCb
        for c = 1:3  % column of aCb
            [eps, k] = perpar([r c 0]);
            aCb(r,c) = kron_delta(r,c)*cosd(theta) - ...
                       eps*lambda(k)*sind(theta) + ...
                       lambda(r)*lambda(c)*(1 - cosd(theta));
        end
    end
end