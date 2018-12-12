%%%%%
% Jordan Mayer
% AAE 532
% 
% f_g_deltaTheta:
%   Get new r, v vector using initial r, v vectors, old and new true
%   anomalies, orbital characteristics, and f and g functions
%
%   Inputs:
%     r_0: original position vector (km)
%     v_0: original velocity vector (km/s)
%     deltaTheta: change in true anomaly (rad)
%     p: semilatus rectum (km)
%     mu: gravitational parameter of central body (km^3/s^2)
%     r_mag: magnitude of new position vector (km)
%
%   Outputs:
%     r: new position vector (km)
%     v: new velocity vector (km/s)
%%%%%

function [r, v] = f_g_deltaTheta(r_0, v_0, deltaTheta, p, mu, r_mag)
    r_0_mag = norm(r_0);  % km
    
    f = 1 - r_mag/p * (1 - cos(deltaTheta));
    g = r_mag * r_0_mag/sqrt(mu*p) * sin(deltaTheta);
    r = f*r_0 + g*v_0;  % km
    
    f_prime = dot(r_0, v_0)/(p*r_0_mag) * (1 - cos(thetaStar)) - ...
              1/r_0 * sqrt(mu/p) * sin(deltaTheta);
    g_prime = 1 - r_0/p * (1 - cos(deltaTheta));
    v = f_prime*r_0 + g_prime*v_0;  % km/s
end