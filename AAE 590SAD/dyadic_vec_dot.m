%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% dyadic_vec_dot:
%   Evaluate dot product between dyadic and vector, assuming expressed in
%   same reference frame.
%
% Inputs:
%   D: dyadic, matrix form
%   v: vector
%
% Outputs:
%   dp: result of D <dot> v
%%%%%

function [dp] = dyadic_vec_dot(D, v)
    dp = zeros(1,3);

    for p = 1:3  % row in D
        for q = 1:3  % column in D
            for e = 1:3  % element in v
                dp(p) = dp(p) + kron_delta(q,e)*D(p,q)*v(e);
            end
        end
    end
end