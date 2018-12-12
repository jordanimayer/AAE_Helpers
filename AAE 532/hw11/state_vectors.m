%%%%%
% Jordan Mayer
% AAE 532
%
% state_vectors:
%   Get position and velocity vectors of an orbiting body on a certain date
%   using orbit characteristics and change in mean anomaly.
%
% Inputs:
%   mu: gravitational parameter of orbited body, m^3/s^2
%   a: semi-major axis of orbit, km
%   e: eccentricity of orbit, km
%   i: inclination of orbit, rad
%   raan: right ascension of ascending node of orbit, rad
%   aop: argument of periapsis of orbit, rad
%   t_0: reference date, JD
%   M_0: reference mean anomaly, rad
%   t_1: date, JD
%
% Outputs:
%   r: position vector in inertial coordinates, km, [x y z]'
%   v: velocity vector in inertial coordinates, km/s, [x y z]'
%%%%%

function [r, v] = state_vectors(mu, a, e, i, raan, aop, t_0, M_0, t_1)
    deltaT = t_1 - t_0  % change in time, days
    deltaT_s = deltaT*24*60*60  % seconds
    n = sqrt(mu/a^3)  % mean motion, rad/s
    deltaM = n*deltaT_s  % change in mean anomaly, rad
    M = M_0 + deltaM  % mean anomaly, rad
    M = zero_to_twoPi(M)  % rad, between 0 and 2pi
    E = E_from_M(e, M)  % eccentric anomaly, rad
    thetaStar = 2*atan(((1+e)/(1-e))^(1/2) * tan(E/2))  % true anomaly, rad
    theta = aop + thetaStar  % rad
    theta = zero_to_twoPi(theta)  % rad, between 0 and 2pi
    
    fprintf('\n');
    
    p = a * (1 - e^2)  % semi-latus rectum, km
    r = p/(1 + e*cos(thetaStar))  % distance, km
    v = sqrt(mu * (2/r - 1/a))  % speed, km/s
    h = sqrt(mu*p)  % angular momentum, km^2/s
    gamma = acos(h/(r*v))  % flight path angle, rad
    if E > pi
        gamma = -gamma  % descending
    end
    v_r = v*sin(gamma)  % radial velocity, km/s
    v_theta = v*cos(gamma)  % tangential velocity, km/s
    
    fprintf('\n');
    
    r_rth = [r 0 0]'  % position vector in r-theta-h coordinates
    v_rth = [v_r v_theta 0]'  % velocity vector in r-theta-h coordinates
    
    fprintf('\n');
    
    r = rthetah_to_xyz(r_rth, i, raan, theta)
      % position vector in x-y-z coordinates
    v = rthetah_to_xyz(v_rth, i, raan, theta)
      % velocity vector in x-y-z coordinates
end
