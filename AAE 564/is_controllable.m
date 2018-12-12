%%%%%
% Jordan Mayer
% AAE 564
%
% is_controllable:
%   Check if (A, B) pair is controllable using PBH test
%
%   Inputs:
%     A, B
%
%   Outputs:
%     controllable: true or false
%%%%%

function [controllable] = is_controllable(A, B)
    controllable = true;
    n = length(A);
    lambda = eig(A);  % eigenvalues of A
    
    for k = 1:length(lambda)
        controllable = (rank([A - lambda(k)*eye(n), B]) == n);
        if (~controllable)
            break;
        end
    end
end