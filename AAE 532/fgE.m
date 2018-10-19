%%%%%
% Jordan Mayer
% AAE 532
% 
% fgE:
%   Get new r, v vector using initial r, v vectors, old and new eccentric
%   anomalies, orbital characteristics, and f and g functions
%
%   Inputs:
%     r_0: original position vector (km)
%     v_0: original velocity vector (km/s)
%     E_0: original eccentric anomaly (rad)
%     t_0_minus_t_p: original time since periapsis (s)
%     E: new eccentric anomaly (rad)
%     a: semimajor axis (km)
%     mu: gravitational parameter of central body (km^3/s^2)
%     e: eccentricity
%
%   Outputs:
%     r: new position vector (km)
%     v: new velocity vector (km/s)
%%%%%

function [r, v] = fgE(r_0, v_0, E_0, t_0_minus_t_p, E, a, mu, e)
  f = 1 - a/norm(r_0) * (1 - cos(E - E_0))
  
  t_minus_t_p = sqrt(a^3 / mu) * (E - e*sin(E));
  t_minus_t_0 = t_minus_t_p - t_0_minus_t_p;
  
  g = t_minus_t_0 - sqrt(a^3 / mu) * ((E - E_0) - sin(E - E_0))
  
  r = f*r_0 + g*v_0;
  
  f_prime = -sqrt(mu * a) / (norm(r) * norm(r_0)) * sin(E - E_0)
  g_prime = 1 - a / norm(r) * (1 - cos(E - E_0))
  
  v = f_prime*r_0 + g_prime*v_0;
end