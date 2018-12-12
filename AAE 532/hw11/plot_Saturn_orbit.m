%%%%%
% Jordan Mayer
% AAE 532
%
% plot_Saturn_orbit:
%   Plot Saturn's orbit in 3-D.
%%%%%

function [] = plot_Saturn_orbit()
    global a_Saturn e_Saturn i_Saturn raan_Saturn aop_Saturn;

    plot_orbit_3D(a_Saturn, e_Saturn, i_Saturn, raan_Saturn, ...
        aop_Saturn, [0 2*pi]);
end
