%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% convert_dyadic_dircos:
%   Convert dyadic from dextral, orthonormal triad vector basis to another
%   using the direction cosine matrix between the two frames.
%
% Inputs:
%   nCb: direction cosine matrix from n-frame to b-frame
%   Ib: dyadic in b-frame, matrix form
%
% Outputs:
%   In: dyadic in n-frame, matrix form
%%%%%

function [In] = convert_dyadic_dircos(nCb, Ib)
    In = zeros(3,3);

    for p = 1:3  % p = row in In
        for q = 1:3  % q = column in In
            for r = 1:3  % r = row in Ib
                for c = 1:3  % c = column in Ib
                    In(p,q) = In(p,q) + nCb(p,r)*Ib(r,c)*nCb(q,c);
                end
            end
        end
    end
end