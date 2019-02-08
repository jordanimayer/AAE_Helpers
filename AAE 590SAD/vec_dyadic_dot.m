%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% vec_dyadic_dot:
%   Evaluate dot product between vector and dyadic, assuming both in same
%   reference frame.
%
% Inputs:
%   v: vector
%   D: dyadic, matrix form
%
% Outputs:
%   dp: result of v <dot> D
%%%%%

function [dp] = vec_dyadic_dot(v, D)
    dp = zeros(1,3);
    
    for p = 1:3
        for q = 1:3
            for e = 1:3
                dp(q) = dp(q) + kron_delta(p,e)*v(e)*D(p,q);
            end
        end
    end
end