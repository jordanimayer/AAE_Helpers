%%%%%
% Jordan Mayer
% AAE 532
%
% plot_Jupiter_orbit:
%   Plot Jupiter's orbit in 3-D.
%%%%%

function [] = plot_Jupiter_orbit()
    global a_Jupiter e_Jupiter i_Jupiter raan_Jupiter aop_Jupiter;

    plot_orbit_3D(a_Jupiter, e_Jupiter, i_Jupiter, raan_Jupiter, ...
        aop_Jupiter, [0 2*pi]);
end
