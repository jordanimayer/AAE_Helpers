%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% perpar:
%   Calculate permutation parameter (also known as Levi-Civita symbol) for
%   three numbers 1, 2, 3 in some order.
%
% Inputs:
%   ijk: numbers (1, 2, or 3), in order: [i j k]
%
% Outputs:
%   eps: permutation parameter for i, j, k
%   k: third number (in case k was unknown)
%%%%%

function [eps, k] = perpar(ijk)
    i = ijk(1);
    j = ijk(2);
    k = ijk(3);
    
    if (k == 0)  % indicates k unknown
        k = perpar_thirdnum(ijk(1:2));
    end
    
    if (i == j || j == k || i == k)
        eps = 0; return;
    elseif (i == 1)
        if (j == 2)
            eps = 1; return;
        else
            eps = -1; return;
        end
    elseif (i == 2)
        if (j == 3)
            eps = 1; return;
        else
            eps = -1; return;
        end
    else
        if (j == 1)
            eps = 1; return;
        else
            eps = -1; return;
        end
    end
end