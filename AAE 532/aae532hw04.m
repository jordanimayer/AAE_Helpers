%%%%%
% Jordan Mayer
% AAE 532
% HW 04
%
% Perform various arithmetic and plot results
%%%%%

%%% Problem 1 %%%

fprintf('\nProblem 1:\n\n');

% Part a: calculate various orbital characteristics and initial
%         spacecraft characteristics

fprintf('Part a:\n\n');

r_earth = 6378.137;  % Earth radius, km (given in table of constants)
mu_earth = 398600.4418;  % Earth gravitational parameter, km^3/s^2
                         % (given in table of constants)
e = 0.4;  % eccentricity (given)
a = 6*r_earth  % semimajor axis, km
E_0 = -90;  % initial eccentric anomaly, deg

p = a * (1-e*e);  % semi-latus rectum, km
r_p = a * (1-e);  % radius at periapsis, km
r_a = a * (1+e);  % radius at apoapsis, km
P = 2*pi * sqrt(a*a*a/mu_earth);  % orbital period, s
[P_hr, P_min, P_s] = hr_min_sec(P);
curlE = -mu_earth/(2*a);  % specific mechanical energy, km^2/s^2

fprintf('p = %.0f km\n', p);
fprintf('r_p = %.0f km\n', r_p);
fprintf('r_a = %.0f km\n', r_a);
fprintf('P = %.0f s = %d hr, %d min, %d s\n', ...
        P, P_hr, P_min, P_s);
fprintf('curlE = %.4f km^2/s^2\n', curlE);

fprintf('\n');

thetaStar_0 = 2*atand(((1+e)/(1-e))^(1/2) * tand(E_0/2));
  % initial true anomaly, deg
r_0 = p/(1+e*cosd(thetaStar_0)); % initial radius, km
v_0 = (2*(curlE + mu_earth/r_0))^(1/2);  % initial speed, km/s
gamma_0 = -asind(a*e/r_0);  % initial flight path angle, deg

fprintf('r_0 = %.0f km\n', r_0);
fprintf('v_0 = %.4f km/s\n', v_0);
fprintf('thetaStar_0 = %.2f deg\n', thetaStar_0);
fprintf('gamma_0 = %.3f deg\n', gamma_0);

% Part b: calculate initial and final conditions

fprintf('\nPart b:\n\n');

r_0_1 = -r_0*sind(abs(gamma_0));  % position vector component in e_hat dir
r_0_2 = -r_0*cosd(abs(gamma_0));  % position vector component in p_hat dir
fprintf('r_0 vector = %.0fe_hat + %.0fp_hat (km)\n', r_0_1, r_0_2);
fprintf('v_0_vector = %.4fe_hat (km/s)\n', v_0);

fprintf('\n');

thetaStar_f = 90;  % final true anomaly, deg
r_f = p;  % final radius, km
v_f = (2*(curlE + mu_earth/r_f))^(1/2);  % final speed, km/s
v_f_1 = v_f*-cosd(gamma_f);  % velocity component in e_hat dir
v_f_2 = v_f*sind(gamma_f);  % velocity component in p_hat dir
E_f = 2*atand(((1-e)/(1+e))^(1/2));  % final eccentric anomaly, deg
gamma_f = 90 - asind(-v_f_1/v_f);  % final flight path angle, deg
fprintf('r_f = %.0f km\n', r_f);
fprintf('v_f = %.4f km/s\n', v_f);
fprintf('E_f = %.3f deg\n', E_f);
fprintf('gamma_f = %.3f deg\n', gamma_f);

fprintf('\n');

fprintf('r_f_vector = %.0fp_hat\n', r_f);
fprintf('v_f_vector = %.4fe_hat + %.4fp_hat\n', v_f_1, v_f_2);

% Part c: calculate times and changes

fprintf('\nPart c:\n\n');

t_0 = sqrt(a*a*a/mu_earth) * (deg2rad(E_0) + e);
  % initial time, s (relative to periapsis)
t_f = sqrt(a*a*a/mu_earth) * (deg2rad(E_f) - e*sind(E_f));
  % final time, s (relative to periapsis)
tof = t_f - t_0;
[t_0_hr, t_0_min, t_0_s] = hr_min_sec(t_0);
[t_f_hr, t_f_min, t_f_s] = hr_min_sec(t_f);
[tof_hr, tof_min, tof_s] = hr_min_sec(tof);
fprintf('t_0 = %.4d s = %d hr, %d min, %d s\n', ...
        t_0, t_0_hr, t_0_min, t_0_s);
fprintf('t_f = %.4d s = %d hr, %d min, %d s\n', ...
        t_f, t_f_hr, t_f_min, t_f_s);
fprintf('tof = %.4d s = %d hr, %d min, %d s\n', ...
        tof, tof_hr, tof_min, tof_s);
fprintf('delta_thetaStar = %.2f deg\n', thetaStar_f - thetaStar_0);
fprintf('delta_E = %.2f deg\n', E_f - E_0);

% Part d: plot entire orbit (to be annotated by hand)

thetaStar = linspace(0,360,1000);  % range of true anomalies, deg
r = p ./ (1 + e*cosd(thetaStar));  % corresponding radii, km
r_1 = r .* cosd(thetaStar);  % radii components in e_hat dir
r_2 = r .* sind(thetaStar);  % radii components in p_hat dir

plot(r_1, r_2); grid on;
xlabel('r_1 (km)'); ylabel('r_2 (km)');
title('Jordan Mayer - HW 04 Problem 1 Part d');
axis([-60000 60000 -60000 60000]);  % scale axes for more intuitive plot


%%% Problem 3 %%%

fprintf('\nProblem 3:\n\n');

% Part a: calculate orbital characteristics

fprintf('Part a:\n\n');

r_mars = 3396.19;  % Mars radius, km (given in table of constants)
mu_mars = 42828.3719012840;  % Mars gravitational parameter, km^3/s^2
                             % (given in table of constants)
v_inf = 2.65;  % incoming speed, km/s
r_p = r_mars + 500;  % radius at periapsis, km

a = -mu_mars/(v_inf*v_inf);  % semimajor axis, km
e = 1 - r_p/a;  % eccentricity
p = a * (1 - e*e);  % semi-latus rectum, km
h = sqrt(mu_mars * p);  % angular momentum, km^2/s
curlE = -mu_mars/(2*a);  % specific mechanical energy, km^2/s^2
delt = 2*asind(1/e);  % flyby angle, deg
thetaStar_inf = acosd(-1/e);  % outgoing true anomaly, deg

fprintf('p = %.0f km\n', p);
fprintf('a = %.1f km\n', a);
fprintf('e = %.4f\n', e);
fprintf('h = %.0f km^2/s\n', h);
fprintf('curlE = %.4f km^2/s^2\n', curlE);
fprintf('delt = %.3f deg\n', delta);
fprintf('thetaStar_inf = %.2f deg\n', thetaStar_inf);

% Part b: calculate semiminor axis and aim point characteristics

fprintf('\nPart b:\n\n');

b = (r_p + abs(a)) * sind(90 - delta/2);  % semiminor axis, km
fprintf('b = %.1f km\n', b);

fprintf('\n');

thetaStar_b = delta/2;  % true anomaly at aim point, deg
r_b = p/(1 + e*cosd(thetaStar_b));  % radius at aim point, km
v_b = (2*(curlE + mu_mars/r_b))^(1/2);  % speed at aim point, km/s
gamma_b = 90 - asind(sqrt(mu_mars*p)/(r_b*v_b));
  % flight path angle at aim point, deg

fprintf('r_b = %.1f km\n', r_b);
fprintf('v_b = %.4f km/s\n', v_b);
fprintf('thetaStar_b = %.3f deg\n', thetaStar_b);
fprintf('gamma_b = %.3f deg\n', gamma_b);

thetaStar = linspace(-thetaStar_inf,thetaStar_inf,1000);
  % range of true anomalies, deg
r = p ./ (1 + e*cosd(thetaStar));  % corresponding radii, km
r_1 = r .* cosd(thetaStar);  % radii components in e_hat dir
r_2 = r .* sind(thetaStar);  % radii components in p_hat dir

% create two plots: one to see asymptotes, one to see closer
% characteristics
plot(r_1, r_2); grid on;
xlabel('r_1 (km)'); ylabel('r_2 (km)');
title('Jordan Mayer - HW 04 Problem 3 Part b');
axis([-2 2 -2 2] * 10^4);
  % scale axes for more intuitive plot
  
% Part d: calculate various speeds related to periapsis

fprintf('\nPart d:\n\n');

v_p = (2 * (curlE + mu_mars/r_p))^(1/2);  % speed at periapsis, km/s
v_c_p = sqrt(mu_mars/r_p);  % speed for circular orbit with r = r_p, km/s
deltaV_p = v_c_p - v_p;  % change in speed to drop into circular orbit with
                         % r = r_p, km/s
fprintf('v_p = %.4f km/s\n', v_p);
fprintf('v_c_p = %.4f km/s\n', v_c_p);
fprintf('deltaV_p = %.4f km/s\n', deltaV_p);