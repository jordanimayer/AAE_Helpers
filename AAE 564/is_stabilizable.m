%%%%%
% Jordan Mayer
% AAE 564
%
% is_stabilizable:
%   Check if (A, B) pair is stabilizable using PBH test
%
%   Inputs:
%     A, B
%
%   Outputs:
%     stabilizable: true or false
%%%%%

function [stabilizable] = is_stabilizable(A, B)
    stabilizable = true;
    n = length(A);
    lambda = eig(A);  % eigenvalues of A
    
    for k = 1:length(lambda)
        if real(lambda(k)) >= 0
            stabilizable = (rank([A - lambda(k)*eye(n), B]) == n);
            if (~stabilizable)
                break;
            end
        end
    end
end