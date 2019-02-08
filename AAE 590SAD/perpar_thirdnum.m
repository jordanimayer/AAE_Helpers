%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% perpar_thirdnum:
%   Calculate third number for permutation parameter, given the first two.
%
% Inputs:
%   ij: first two numbers (1, 2, or 3) in order: [i j]
%
% Outputs:
%   k: third number
%%%%%

function [k] = perpar_thirdnum(ij)
    i = ij(1);
    j = ij(2);
    
    if (i == j)
        k = i; return;
    elseif (i == 1)
        if (j == 2)
            k = 3; return;
        else
            k = 2; return;
        end
    elseif (i == 2)
        if (j == 3)
            k = 1; return;
        else
            k = 3; return;
        end
    else
        if (j == 1)
            k = 2; return;
        else
            k = 1; return;
        end
    end
end