%%%%%
% Jordan Mayer
% AAE 532
%
% plotOrbit2D:
%   Create a 2-D plot of an orbit, where x is radius in direction
%   periapsis and y is radius in direction of semi-latus rectum
%
%   inputs:
%       e: eccentricity
%       p: semi-latus rectum (km)
%       thetaStar_range: range of true anomalies to plot (deg)
%                        [thetaStar_min, thetaStar_max]
%       ttl: string with title of plot
%%%%%

function [] = plotOrbit2D(e,p, thetaStar_range, ttl)
    thetaStar = linspace(thetaStar_range(1), thetaStar_range(2) ,1000);
      % range of true anomalies, deg
    r = p ./ (1 + e*cosd(thetaStar));  % corresponding radii, km
    r_1 = r .* cosd(thetaStar);  % radii components in e_hat dir
    r_2 = r .* sind(thetaStar);  % radii components in p_hat dir

    plot(r_1, r_2); grid on;
    xlabel('r_1 (km)'); ylabel('r_2 (km)');
    title(ttl);
end
