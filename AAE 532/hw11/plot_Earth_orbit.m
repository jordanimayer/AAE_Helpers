%%%%%
% Jordan Mayer
% AAE 532
%
% plot_Earth_orbit:
%   Plot Earth's orbit in 3-D.
%%%%%

function [] = plot_Earth_orbit()
    global a_Earth e_Earth i_Earth raan_Earth aop_Earth;

    plot_orbit_3D(a_Earth, e_Earth, i_Earth, raan_Earth, aop_Earth, ...
        [0 2*pi]);
end
