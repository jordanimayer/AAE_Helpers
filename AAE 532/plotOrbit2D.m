%%%%%
% Jordan Mayer
% AAE 532
%
% plotOrbit2D:
%   Create a 2-D plot of an orbit, where x is radius in direction
%   original periapsis and y is radius in direction of original
%   semi-latus rectum
%
%   inputs:
%       e: eccentricity
%       p: semi-latus rectum (km)
%       thetaStar_range: range of true anomalies to plot (deg)
%                        [thetaStar_min, thetaStar_max]
%       aop: argument of periapsis (deg), 0 for original orbit
%       ttl: string with title of plot
%%%%%

function [] = plotOrbit2D(e, p, thetaStar_range, aop, ttl)
    thetaStar = linspace(thetaStar_range(1), thetaStar_range(2), 1000);
      % range of true anomalies, deg
    aop = ones(1,1000) * aop;
    r = p ./ (1 + e*cosd(thetaStar));  % corresponding radii, km
    r_1 = r .* cosd(thetaStar + aop);  % radii components in e_hat dir
    r_2 = r .* sind(thetaStar + aop);  % radii components in p_hat dir

    plot(r_1, r_2); grid on;
    xlabel('r_1 (km)'); ylabel('r_2 (km)');
    title(ttl);
end
