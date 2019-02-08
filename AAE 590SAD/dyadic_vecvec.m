%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% dyadic_vecvec:
%   Form dyadic from two vectors (expressed in the same frame).
%
% Inputs:
%   v: first vector
%   u: second vector
%
% Outputs:
%   D: resulting dyadic
%%%%%

function [D] = dyadic_vecvec(v, u)
    D = zeros(3,3);
    
    for e1 = 1:3  % element of v1
        for e2 = 1:3  % element of v2
            D(e1,e2) = D(e1,e2) + v(e1)*u(e2);
        end
    end
end