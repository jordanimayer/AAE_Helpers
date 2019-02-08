%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% dyadic_vec_cross:
%   Perform cross product between dyadic and vector, assuming they are
%   expressed in the same frame.
%
% Inputs:
%   D: dyadic, matrix form
%   v: vector
%
% Outputs:
%   C: result of D <cross> v
%%%%%

function [C] = dyadic_vec_cross(D, v)
    C = zeros(3,3);
    
    for p = 1:3  % row of D
        for q = 1:3  % column of D
            for e = 1:3  % element of v
                [eps, k] = perpar([q e 0]);
                C(p,k) = C(p,k) + eps*D(p,q)*v(e);
            end
        end
    end
end