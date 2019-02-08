%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% dircos_b212:
%   Calculate direction cosine matrix from frame A to frame B, where frame
%   B is initially frame A and then is rotated in a Body-two 2-1-2
%   sequence.
%
% Inputs:
%   phi: vector of rotation angles [phi1 phi2 phi3], in 2-1-2 sequence, deg
%
% Outputs:
%   C: direction cosine matrix
%%%%%

function [C] = dircos_b212(phi)
    s1 = sind(phi(1));
    s2 = sind(phi(2));
    s3 = sind(phi(3));
    c1 = cosd(phi(1));
    c2 = cosd(phi(2));
    c3 = cosd(phi(3));

    C = zeros(3,3);
    
    C(1,1) = -s1*c2*s3 + c3*c1;
    C(1,2) = s1*s2;
    C(1,3) = s1*c2*c3 + s3*c1;
    C(2,1) = s2*s3;
    C(2,2) = c2;
    C(2,3) = -s2*c3;
    C(3,1) = -c1*c2*s3 - c3*s1;
    C(3,2) = c1*s2;
    C(3,3) = c1*c2*c3 - s3*s1;
end