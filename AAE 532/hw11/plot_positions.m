%%%%%
% Jordan Mayer
% AAE 532
%
% plot_positions:
%   Plot specific orbital positions in 3-D.
%
% Inputs:
%   r_list: array of vertical x-y-z position vectors, km, [r1 r2 r3 ...]
%%%%%

function [] = plot_positions(r_list)
    for k=1:size(r_list,2)
        r = r_list(:,k);
        plot3(r(1), r(2), r(3), '*');
    end
end
