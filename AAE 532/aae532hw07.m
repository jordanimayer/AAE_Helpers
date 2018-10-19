%%%%%
% Jordan Mayer
% AAE 532
% HW 07
%%%%%

clear all; close all; format long; format compact;

%%% Problem 2 %%%

fprintf('\nProblem 2:\n');

fprintf('\nPart c:\n\n');

thetaStar1_minus = 150;  % true anomaly in original orbit at maneuver
                         % point, deg
raan_minus = 60;  % right ascension of ascending node in original orbit,
                  % deg
omega_minus = 90;  % argument of periapsis in original orbit, deg
i_minus = 30;  % inclination of original orbit, deg
r1_rth0 = [62702.9 0 0]';
  % position vector at maneuver point in original r-theta-h coordinates, km
v1_plus_rth0 = [1.82621 2.09175 -0.5]';
  % velocity vector at maneuver point immediately after maneuver in
  % original r-theta-h coordinates, km/s

DCM_0 = DCM_rthetah_xyz(raan_minus, i_minus, omega_minus, thetaStar1_minus)

r_hat_0 = DCM_0(:,1)
  % radial unit vector of original orbit
theta_hat_0 = DCM_0(:,2)
  % tangential unit vector of original orbit
h_hat_0 = DCM_0(:,3)
  % angular momentum unit vector of original orbit

r1_xyz = DCM_0 * r1_rth0
  % position vector at maneuver point in inertial coordinates, km
v1_plus_xyz = DCM_0 * v1_plus_rth0
  % velocity vector at maneuver point immediately after maneuver in
  % inertial coordinates, km/s
gamma1_plus = acosd(dot(r1_xyz, v1_plus_xyz)/...
                    (norm(r1_xyz) * norm(v1_plus_xyz)))
  % flight path angle at maneuver point immediately after maneuver, deg
  % result confirms previous work
  
fprintf('\nPart d:\n\n');

h_hat_N = cross(r1_xyz, v1_plus_xyz)/norm(cross(r1_xyz, v1_plus_xyz))
  % angular momentum unit vector of new orbit
r_hat_N = DCM_0(:,1)
  % radial unit vector of new orbit
theta_hat_N = cross(h_hat_N, r_hat_N)
  % local horizon unit vector of new orbit
  
DCM_N = [r_hat_N theta_hat_N h_hat_N]

r1_rthN = inv(DCM_N) * r1_xyz
  % position vector at maneuver point in new r-theta-h coordinates, km
v1_plus_rthN = inv(DCM_N) * v1_plus_xyz
  % velocity vector at maneuver point immediately after maneuver in
  % new r-theta-h coordinates, km/s

fprintf('\nPart f:\n\n');

r_0 = r1_xyz;
  % inertial position vector immediately following maneuver, km
v_0 = v1_plus_xyz;
  % inertial velocity vector immediately following maneuver, km/s
E_0 = deg2rad(45.5179);
  % eccentric anomaly immediately following maneuver, rad
t_0_minus_t_p = 21664.1;
  % time since periapsis immediately following maneuver, s
E = deg2rad(90);  % eccentric anomaly at new position, rad
a = 135270;  % semimajor axis of new orbit, km
mu = 324859;  % gravitational parameter of Venus, km^3/s^2
e = 0.765624;  % eccentricity of new orbit
  
[r, v] = fgE(r_0, v_0, E_0, t_0_minus_t_p, E, a, mu, e)

%%% Problem 3 %%%

fprintf('\nProblem 3:\n');

% Part a

fprintf('\nPart a:\n');

% Part i

fprintf('\nPart i:\n\n');

R_earth = 6378.14;  % Earth radius, km
mu_earth = 398600;  % Earth gravitational parameter, km^3/s^2
r1 = 1.04 * R_earth;  % initial orbit radius, km
r2 = 6.6 * R_earth;  % final orbit radius, km

% first maneuver

fprintf('\nFirst maneuver:\n\n');

v1_minus = sqrt(mu_earth/r1)
  % velocity magnitude at first burn point before maneuver, km/s
curlE_T = -mu_earth/(r1 + r2)
  % specific mechanical energy of transfer orbit, km^2/s^2
v1_plus = sqrt(2 * (curlE_T + mu_earth/r1))
  % velocity magnitude at first burn point after maneuver, km/s
deltaV_1 = v1_plus - v1_minus
  % deltaV magnitude for first maneuver, km/s

% second maneuver

fprintf('\nSecond maneuver:\n\n');

v2_minus = sqrt(2 * (curlE_T + mu_earth/r2))
  % velocity magnitude at second burn point before maneuver, km/s
v2_plus = sqrt(mu_earth/r2)
  % velocity magnitude at second burn point after maneuver, km/s
deltaV_2 = v2_plus - v2_minus
  % deltaV magnitude for second maneuver, km/s

% transfer characteristics

fprintf('\nTransfer characteristics:\n\n');

deltaV_Tot = deltaV_1 + deltaV_2
  % total transfer deltaV magnitude, km/s
a_T = (r1 + r2)/2
  % semimajor axis of transfer orbit, km
IP_T = 2*pi*sqrt(a_T^3/mu_earth)
  % period of transfer orbit, s
TOF = 0.5*IP_T/3600
  % time of flight of transfer orbit, hr
  
% for plotting in part b
e_T_i = r2/a_T - 1;
p_T_i = r2 * (1 - e_T_i);

% for part c
TOF_i = TOF;
  
% Part ii

fprintf('\nPart ii:\n');

% first maneuver

fprintf('\nFirst maneuver:\n\n');

thetaStar_T_1 = 60;
  % true anomaly of transfer orbit at first burn point, deg
e_T = (r2 - r1)/(r1*cosd(thetaStar_T_1) + r2)
  % eccentricity of transfer orbit
a_T = r2/(1 + e_T)
  % semimajor axis of transfer orbit, km
curlE_T = -mu_earth/(2*a_T)
  % specific mechanical energy of transfer orbit, km^2/s^2
v1_plus = sqrt(2 * (curlE_T + mu_earth/r1))
  % velocity magnitude at first burn point after maneuver, km/s
p_T = r2 * (1 - e_T)
  % semi-latus rectum of transfer orbit, km
h_T = sqrt(mu_earth*p_T)
  % angular momentum of transfer orbit, km^2/s
gamma1_plus_1 = acosd(h_T/(r1 * v1_plus))
gamma1_plus_2 = -gamma1_plus_1
  % deg, two options to account for quadrant ambiguity
gamma1_plus = gamma1_plus_1;
  % flight path angle at first burn point after maneuver, deg
  % choose positive option b/c ascending
deltaV_1 = sqrt(v1_minus^2 + v1_plus^2 - ...
                2*v1_minus*v1_plus*cosd(gamma1_plus))
  % deltaV for first maneuver, km/s
beta1_1 = asind(v1_plus/deltaV_1 * sind(gamma1_plus))
beta1_2 = 180 - beta1_1
  % deg, two options to account for quadrant ambiguity
eta1_1 = asind(v1_minus/deltaV_1 * sind(gamma1_plus))
eta1_2 = 180 - eta1_1
  % deg, two options to account for quadrant ambiguity
triangle11 = gamma1_plus + beta1_1 + eta1_1
triangle12 = gamma1_plus + beta1_1 + eta1_2
triangle21 = gamma1_plus + beta1_2 + eta1_1
triangle22 = gamma1_plus + beta1_2 + eta1_2
  % correct combination should sum to 180 deg
alpha1 = 180 - beta1_2
  % turn angle of first maneuver, deg
  
% second maneuver

fprintf('\nSecond maneuver:\n\n');

v2_minus = sqrt(2 * (curlE_T + mu_earth/r2))
  % velocity magnitude at second burn point before maneuver
deltaV_2 = v2_plus - v2_minus
  % deltaV magnitude of second maneuver
  
% transfer characteristics

fprintf('\nTransfer characteristics:\n\n');

IP_T = 2*pi*sqrt(a_T^3 / mu_earth)
  % transfer orbit period, s
TOF = 0.5 * IP_T / 3600
  % time of flight for transfer, hr
deltaV_Tot = deltaV_1 + deltaV_2
  % total deltaV magnitude for transfer, km/s
  
% for plotting in part b
e_T_ii = e_T;
p_T_ii = p_T;
  
% Part b

fprintf('\nPart b:\n\n');

% for transfer ii
deltaV_1_r = deltaV_1 * sind(alpha1)
  % radial deltaV for first maneuver, km/s
deltaV_1_theta = deltaV_1 * cosd(alpha1)
  % tangential deltaV for first maneuver, km/s
  
% create plots

% transfer i
figure();
hold on;
ttl = 'Jordan Mayer - HW 07 Problem 3b, Transfer i (Hohmann)';
plotOrbit2D(0, r1, [0 360], '');
plotOrbit2D(e_T_i, p_T_i, [0 180], '');
plotOrbit2D(0, r2, [0 360], ttl);
legend('Original Orbit', 'Transfer Orbit', 'Final Orbit');

% transfer ii
figure();
hold on;
ttl = ['Jordan Mayer - HW 07 Problem 3b, Transfer ii' ...
       '(non-tangential departure)'];
plotOrbit2D(0, r1, [0 360], '');
plotOrbit2D(e_T_ii, p_T_ii, [thetaStar_T_1 180], '');
plotOrbit2D(0, r2, [0 360], ttl);
legend('Original Orbit', 'Transfer Orbit', 'Final Orbit');

% Part c

fprintf('\nPart c:\n\n');
TOF_i = TOF_i * 3600
  % TOF of Hohmann transfer, s
vc_2 = v2_plus;
  % circular velocity at radius 2 (velocity of rendezvous craft), km/s
n_rc = rad2deg(vc_2/r2)
  % angular velocity of rendezvous craft, deg/s
phi = 180 - n_rc * TOF_i
  % phase angle required at LEO departure, deg
  
%%% Problem 4 %%%

fprintf('\nProblem 4:\n');

% Part a

fprintf('\nPart a:\n');

mu_sun = 1.32712 * 10^11;  % Sun gravitational parameter, km^3/s^2
r1 = 149650000;  % initial orbit radius, km
r2 = 5868130000;  % final orbit radius, km

% first maneuver

fprintf('\nFirst maneuver:\n\n');

v1_minus = sqrt(mu_sun/r1)
  % velocity magnitude at first burn point before maneuver, km/s
curlE_T = -mu_sun/(r1 + r2)
  % specific mechanical energy of transfer orbit, km^2/s^2
v1_plus = sqrt(2 * (curlE_T + mu_sun/r1))
  % velocity magnitude at first burn point after maneuver, km/s
deltaV_1 = v1_plus - v1_minus
  % deltaV magnitude for first maneuver, km/s
  
% second maneuver

fprintf('\nSecond maneuver:\n\n');

v2_plus = sqrt(mu_sun/r2)
  % velocity magnitude at second burn point after maneuver, km/s
v2_minus = sqrt(2 * (curlE_T + mu_sun/r2))
  % velocity magnitude at second burn point before maneuver, km/s
deltaV_2 = v2_plus - v2_minus
  % deltaV magnitude for second maneuver, km/s
  
% transfer characteristics

fprintf('\nTransfer characteristics:\n\n');

deltaV_Tot = deltaV_1 + deltaV_2
  % total deltaV magnitude for transfer, km/s
a_T = (r1 + r2)/2
  % semimajor axis of transfer orbit, km
IP_T = 2*pi*sqrt(a_T^3 / mu_sun)
  % period of transfer orbit, s
TOF = 0.5 * IP_T / 3600 / 24 / 365
  % time of flight of transfer, yr
  
% for part c
TOF_a = TOF * 3600 * 24 * 365;
  % time of flight of transfer, s
v_c_2 = v2_plus;
r_c_2 = r2;
  
% Part b

fprintf('\nPart b:\n\n');

a_2 = 5868130000;  % semimajor axis of final orbit, km
e_2 = 0.247576;  % eccentricity of final orbit
r2 = a_2 * (1 - e_2)
  % perihelion radius of final orbit

% first maneuver

fprintf('\nFirst maneuver:\n\n');

curlE_T = -mu_sun / (r1 + r2)
  % specific mechanical energy of transfer orbit, km^2/s^2
v1_plus = sqrt(2 * (curlE_T + mu_sun/r1))
  % velocity magnitude at first burn point after maneuver, km/s
deltaV_1 = v1_plus - v1_minus
  % deltaV magnitude for first maneuver, km/s

% second maneuver

fprintf('\nSecond maneuver:\n\n');

v2_minus = sqrt(2 * (curlE_T + mu_sun/r2))
  % velocity magnitude at second burn point before maneuver, km/s
curlE_2 = -mu_sun/(2 * a_2)
  % specific mechanical energy of final orbit, km^2/s^2
v2_plus = sqrt(2 * (curlE_2 + mu_sun/r2))
  % velocity magnitude at second burn point after maneuver, km/s
deltaV_2 = v2_plus - v2_minus
  % deltaV magnitude for second maneuver, km/s
  
% transfer characteristics

fprintf('\nTransfer characteristics:\n\n');

deltaV_Tot = deltaV_1 + deltaV_2
  % total deltaV magnitude for transfer, km/s
a_T = (r1 + r2)/2
  % semimajor axis of transfer orbit, km
IP_T = 2*pi*sqrt(a_T^3/mu_sun)
  % period of transfer orbit, s
TOF = 0.5 * IP_T / 3600 / 24 / 365
  % time of flight of transfer, yr

% Part c

fprintf('\nPart c:\n\n');

n_pluto = rad2deg(v_c_2 / r_c_2)
  % angular velocity of pluto in circular orbit, deg/s
phi = 180 - n_pluto * TOF_a
  % phase angle required at departure for rendezvous with Pluto