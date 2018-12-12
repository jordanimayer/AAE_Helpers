%%%%%
% Jordan Mayer
% AAE 532
%
% plot_Neptune_orbit:
%   Plot Neptune's orbit in 3-D.
%%%%%

function [] = plot_Neptune_orbit()
    global a_Neptune e_Neptune i_Neptune raan_Neptune aop_Neptune;

    plot_orbit_3D(a_Neptune, e_Neptune, i_Neptune, raan_Neptune, ...
        aop_Neptune, [0 2*pi]);
end
