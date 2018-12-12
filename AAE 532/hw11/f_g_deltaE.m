%%%%%
% Jordan Mayer
% AAE 532
% 
% f_g_deltaE:
%   Get new r, v vector using initial r, v vectors, old and new eccentric
%   anomalies, orbital characteristics, and f and g functions
%
%   Inputs:
%     r_0: original position vector (km)
%     v_0: original velocity vector (km/s)
%     deltaE: change in eccentric anomaly (rad)
%     deltaT: change in time (s)
%     a: semimajor axis (km)
%     mu: gravitational parameter of central body (km^3/s^2)
%
%   Outputs:
%     r: new position vector (km)
%     v: new velocity vector (km/s)
%%%%%

function [r, v] = f_g_deltaE(r_0, v_0, deltaE, deltaT, a, mu)
    r_0_mag = norm(r_0);  % km
    
    f = 1 - a/r_0_mag * (1 - cos(deltaE));
    g = deltaT - sqrt(a^3 / mu) * (deltaE - sin(deltaE));
    r = f*r_0 + g*v_0;  % km
    r_mag = norm(r);  % km
    
    f_prime = -sqrt(mu*a)/(r_mag*r_0_mag) * sin(deltaE);
    g_prime = 1 - a/r_mag * (1 - cos(deltaE));
    v = f_prime*r_0 + g_prime*v_0;  % km/s
end