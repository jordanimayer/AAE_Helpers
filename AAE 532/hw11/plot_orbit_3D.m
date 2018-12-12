%%%%%
% Jordan Mayer
% AAE 532
%
% plot_orbit_3D:
%   Create 3-dimensional plot of orbit.
%
% Inputs:
%   a: semimajor axis, km
%   e: eccentricity
%   i: inclination, rad
%   raan: right ascension of ascending node, rad
%   aop: argument of periapsis, rad
%   thetatStar_range: range of true anomalies, rad, [min max]
%%%%%

function [] = plot_orbit_3D(a, e, i, raan, aop, thetaStar_range)
    p = a * (1 - e^2);  % semi-latus rectum, km
    x = zeros(1,1000);  % x-position, km
    y = zeros(1,1000);  % y-position, km
    z = zeros(1,1000);  % z-position, km
    thetaStars = linspace(thetaStar_range(1), thetaStar_range(2), 1000);
    
    for k=1:1000
        thetaStar = thetaStars(k);
        
        r_mag = p/(1 + e*cos(thetaStar));  % distance, km
        r_rth = [r_mag 0 0]';  % position vector in r-theta-h coords
        theta = aop + thetaStar;  % rad
        r = rthetah_to_xyz(r_rth, i, raan, theta);
        
        x(k) = r(1);
        y(k) = r(2);
        z(k) = r(3);
    end
    
    max_val = max([x y z]) * 1.05;
    min_val = min([x y z]) * 1.05;
    
    plot3(x, y, z);
    grid on;
    xlim([min_val max_val]);
    ylim([min_val max_val]);
    zlim([min_val max_val]);
    axis square;
    xlabel('r_x (km)');
    ylabel('r_y (km)');
    zlabel('r_z (km)');
end
