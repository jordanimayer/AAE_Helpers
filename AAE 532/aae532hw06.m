%%%%%
% Jordan Mayer
% AAE 532
% HW 06
%
% Create plots of orbits and calculate DCM matrix for 3-D problem.
%%%%%

clear all; close all;

%%% Problem 1 %%%

% Part aiii

e_old = 0.75;  % eccentricity of original orbit (given)
e_new = 0;  % eccentricity of new orbit (circular)
r_1 = 49111.7;  % radius at maneuver point, km (see handwritten work)
p_old = r_1 * (1 - e_old);  % semi-latus rectum of original orbit, km
                            % (maneuver point at apoapsis)
p_new = r_1;  % semi-latus rectum of new orbit, km (circular)

ttl = 'Jordan Mayer - HW 06, Problem 1a, old and new orbits (full)';
figure();
plotOrbit2D(e_old, p_old, [0 360], ttl);
hold on;
plotOrbit2D(e_new, p_new, [0 360], ttl);
axis square;

ttl = 'Jordan Mayer - HW 06, Problem 1a, old and new orbits (maneuver)';
figure();
plotOrbit2D(e_old, p_old, [175 195], ttl);
hold on;
plotOrbit2D(e_new, p_new, [175 195], ttl);
axis square;

% Part biv

r_1 = 42095.7;  % radius at maneuver point, km (see handwritten work)
p_new = r_1;  % semi-latus rectum of new orbit, km (circular)

ttl = 'Jordan Mayer - HW 06, Problem 1b, old and new orbits';
figure();
plotOrbit2D(e_old, p_old, [0 360], ttl);
hold on;
plotOrbit2D(e_new, p_new, [0 360], ttl);
axis square;

ttl = 'Jordan Mayer - HW 06, Problem 1b, old and new orbits (maneuver)';
figure();
plotOrbit2D(e_old, p_old, [-170 -150], ttl);
hold on;
plotOrbit2D(e_new, p_new, [-170 -150], ttl);
axis square;

%%% Problem 2 %%%

% Part d

fprintf('\nProblem 2d:\n\n');

% Get new eccentric anomaly, period, time from periapsis, periapsis radius
e_new = 0.545110;
a_new = 17023.0;  % km
mu_earth = 398600;  % km^3/s^2
thetaStar_new = 124.315;  % deg
E_new = 2*atand(((1 - e_new)/(1 + e_new))^(1/2) * tand(thetaStar_new/2))
  % deg
period_new = 2*pi*sqrt(a_new^3 / mu_earth)  % seconds
time_new = sqrt(a_new^3 / mu_earth) * (deg2rad(E_new) - e_new*sind(E_new))
  % seconds
rp_new = a_new * (1 - e_new)  % m

%%% Problem 3 %%%

% Part a

fprintf('\nProblem 3a:\n\n');

e = 0.5;
M_1 = 4.37652;
E_1 = 123;  % indicates E not found
for E = linspace(0,2*pi,10^6)
    M = E - e*sin(E);
    if abs(M_1 - M) < 0.00001
        E_1 = E
        break;
    end
end

% Part b

fprintf('\nProblem 3b:\n\n');

i = 30;  % inclination, deg
raan = 45;  % right ascension of ascending node, deg
omega = -60;  % argument of periapsis, deg
thetaStar_1 = 209.550;  % true anomaly, deg

r_1_rth = [22539.6 0 0]';  % r_1 in r-theta-h coordinates, km
v_1_rth = [-0.433153 1.04428 0]';  % v_1 in r-theta-h coordinates, km/s
DCM = DCM_rthetah_xyz(raan, i, omega, thetaStar_1);
  % transformation matrix from r-theta-h coordinates to x-y-z (inertial)
  % coordinates
r_1_xyz = DCM * r_1_rth
v_1_xyz = DCM * v_1_rth

% Part c

fprintf('\nProblem 3c:\n\n');

deltaV_xyz = [0.1 -0.2 0.25]';  % deltaV in r-theta-h coordinates, km/s
deltaV_rth = DCM\deltaV_xyz
deltaV_rt = deltaV_rth(1:2)
magDeltaV_rt = norm(deltaV_rt)
magDeltaV = norm(deltaV_rth)
beta = acosd(magDeltaV_rt / magDeltaV)  % deg
magDeltaV_r = deltaV_rth(1)
phi = acosd(magDeltaV_r/magDeltaV_rt)  % deg

v1CrossDeltaV = cross(v_1_rth, [deltaV_rt; 0])
magV1CrossDeltaV = norm(v1CrossDeltaV);
v_1 = norm(v_1_rth);
alpha = asind(magV1CrossDeltaV / (v_1 * magDeltaV_rt))  % deg

% Part d

fprintf('\nProblem 3d:\n\n');

p = 12735.7;  % km
e = 0.5;
a = 5*3396.19;  % km
mu_mars = 42828.4;  % km^3 / s^2
h = 23537.7;  % km^2/s
thetaStar_1 = 240;  % deg

DCM = DCM_rthetah_xyz(raan, i, omega, thetaStar_1);

r_1_mag = p/(1 + e*cosd(thetaStar_1))  % km
v_1_mag = (2 * (mu_mars/r_1_mag - mu_mars/(2*a)))^(1/2)  % km/s
beta_1 = 180 - asind(h/(r_1_mag * v_1_mag))  % deg
gamma_1 = 90 -  beta_1  % deg
v_1 = [v_1_mag * cosd(beta_1); v_1_mag * sind(beta_1); 0]  % km/s
v_1 = DCM * v_1  % in xyz coords
deltaV = [0.1 -0.2 0.25]';  % km/s
v_1_N = v_1 + deltaV  % km/s
v_1_N_mag = norm(v_1_N)
fancE = v_1_N_mag^2 / 2 - mu_mars/r_1_mag  % km^2 / s^2
a_N = -mu_mars / (2*fancE)  % km
r_1 = [r_1_mag 0 0]';
r_1 = DCM * r_1
p_N = (norm(cross(r_1, v_1_N)))^2 / mu_mars  % km
e_N = sqrt(1 - p_N/a_N)