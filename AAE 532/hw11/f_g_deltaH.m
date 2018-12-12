%%%%%
% Jordan Mayer
% AAE 532
% 
% f_g_deltaH:
%   Get new r, v vector using initial r, v vectors, old and new hyperbolic
%   anomalies, orbital characteristics, and f and g functions
%
%   Inputs:
%     r_0: original position vector (km)
%     v_0: original velocity vector (km/s)
%     deltaH: change in eccentric anomaly (rad)
%     deltaT: change in time (s)
%     abs_a: magnitude of semimajor axis (km)
%     mu: gravitational parameter of central body (km^3/s^2)
%
%   Outputs:
%     r: new position vector (km)
%     v: new velocity vector (km/s)
%%%%%

function [r, v] = f_g_deltaH(r_0, v_0, deltaH, deltaT, abs_a, mu)
    r_0_mag = norm(r_0);  % km
    
    f = 1 - abs_a/r_0_mag * (cosh(deltaH) - 1);
    g = deltaT - sqrt(abs_a^3 / mu) * (sinh(deltaH) - deltaH);
    r = f*r_0 + g_v_0;  % km
    r_mag = norm(r);  % km
    
    f_prime = sqrt(mu*abs_a)/(r_mag*r_0_mag) * sinh(deltaH);
    g_prime = 1 - abs_a/r_mag * (cosh(deltaH) - 1);
    v = f_prime*r_0 + g_prime*v_0;  % km/s
end