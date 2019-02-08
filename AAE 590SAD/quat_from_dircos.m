%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% quat_from_dircos:
%   Calculate quaternion between from reference frame A to frame B using
%   direction cosine from A to B. Note that vector component will be
%   valid in both frames.
%
% Inputs:
%   C: direction cosine matrix from frame A to frame B
%
% Outputs:
%   q: quaternion from frame A to frame B, [vector scalar]
%%%%%

function [q] = quat_from_dircos(C)
    q = zeros(1,4);
    
    q(4) = 1/2 * (1 + C(1,1) + C(2,2) + C(3,3))^(1/2);
    q(1) = (C(3,2) - C(2,3))/(4*q(4));
    q(2) = (C(1,3) - C(3,1))/(4*q(4));
    q(3) = (C(2,1) - C(1,2))/(4*q(4));
end