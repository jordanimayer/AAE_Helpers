%%%%%
% Jordan Mayer
% AAE 532
% HW 11
%
% Perform various arithmetic and plotting to perform orbital
% analysis for Pioneer 11 mission.
%%%%%

clear all; close all; format compact; format long; rehash;

%%% Preliminary Setup %%%

global mu_Sun ...
       R_Earth mu_Earth a_Earth e_Earth i_Earth raan_Earth aop_Earth ...
       M0_Earth ...
       R_Jupiter mu_Jupiter a_Jupiter e_Jupiter i_Jupiter raan_Jupiter ...
       aop_Jupiter M0_Jupiter ...
       R_Saturn mu_Saturn a_Saturn e_Saturn i_Saturn raan_Saturn ...
       aop_Saturn M0_Saturn ...
       a_Neptune;
load_constants_planets_moon;  % R, mu of planets
load_planets_11_06_1991;  % a, e, i, raan, aop, M0 of planets on 11/6/1991
fig_position = [200 200 1000 800];  % position and size for figures

% View all planetary orbits
fig_allOrbits = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_Jupiter_orbit;
plot_Saturn_orbit;
plot_Neptune_orbit;
legend({'Sun', 'Earth', 'Jupiter', 'Saturn', 'Neptune'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Planet Orbits');
mkdir('../img/allOrbits');
save_all_views(fig_allOrbits, '../img/allOrbits/allOrbits');

%%%%%%%% Part a: Transfer Orbits %%%%%%%%

fprintf('\n--------Part a: Transfer Orbits--------\n');

%%%%% Part i %%%%%

fprintf('\n-----Part i-----\n\n');

% times calculated using JPL JD calculator
t_1 = 2441779  % JD
t_2 = 2442385  % JD
t_3 = 2444118  % JD
t_4 = 2447946  % JD

%%%%% Part ii %%%%%

fprintf('\n-----Part ii-----\n');

t_0 = 2448566.5;  % reference time (for reference mean anomaly), JD

%%% Earth launch %%%
fprintf('\n---Earth launch---\n\n');

% get Earth state vectors
[r_Earth, v_Earth] = state_vectors(mu_Sun, a_Earth, e_Earth, i_Earth, ...
    raan_Earth, aop_Earth, t_0, M0_Earth, t_1);
r_Earth_1 = r_Earth;  % for later
v_Earth_1 = v_Earth;  % for later

% plot Earth position
fig_earthPos1 = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_positions(r_Earth);
legend({'Sun', 'Earth orbit', 'Earth position'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Earth Position at Launch');
mkdir('../img/earthPos1');
save_all_views(fig_earthPos1, '../img/earthPos1/earthPos1');
                               
% plot all planet positions
[r_Jupiter, dummy] = state_vectors_quiet(mu_Sun, a_Jupiter, e_Jupiter, ...
    i_Jupiter, raan_Jupiter, aop_Jupiter, t_0, M0_Jupiter, t_1);
[r_Saturn, dummy] = state_vectors_quiet(mu_Sun, a_Saturn, e_Saturn, ...
    i_Saturn, raan_Saturn, aop_Saturn, t_0, M0_Saturn, t_1);
fig_allPos1 = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_Jupiter_orbit;
plot_Saturn_orbit;
plot_positions([r_Earth, r_Jupiter, r_Saturn]);
legend({'Sun', 'Earth orbit', 'Jupiter orbit', 'Saturn orbit', ...
       'Earth position', 'Jupiter position', 'Saturn position'}, ...
       'FontSize', 14);
title('Jordan Mayer, HW 11 - Planet Positions at Launch');
mkdir('../img/allPos1');
save_all_views(fig_allPos1, '../img/allPos1/allPos1');

%%% Jupiter encounter %%%
fprintf('\n---Jupiter encounter---\n\n');

% get Jupiter state vectors
[r_Jupiter, v_Jupiter] = state_vectors(mu_Sun, a_Jupiter, e_Jupiter, ...
    i_Jupiter, raan_Jupiter, aop_Jupiter, t_0, M0_Jupiter, t_2);
r_Jupiter_2 = r_Jupiter;  % for later
v_Jupiter_2 = v_Jupiter;  % for later

% plot Jupiter position
fig_jupiterPos2 = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Jupiter_orbit;
plot_positions(r_Jupiter);
legend({'Sun', 'Jupiter orbit', 'Jupiter position'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Jupiter Position at Jupiter Encounter');
mkdir('../img/jupiterPos2');
save_all_views(fig_jupiterPos2, '../img/jupiterPos2/jupiterPos2');

% plot all planet positions
[r_Earth, dummy] = state_vectors_quiet(mu_Sun, a_Earth, e_Earth, ...
    i_Earth, raan_Earth, aop_Earth, t_0, M0_Earth, t_2);
[r_Saturn, dummy] = state_vectors_quiet(mu_Sun, a_Saturn, e_Saturn, ...
    i_Saturn, raan_Saturn, aop_Saturn, t_0, M0_Saturn, t_2);
fig_allPos2 = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_Jupiter_orbit;
plot_Saturn_orbit;
plot_positions([r_Earth, r_Jupiter, r_Saturn]);
legend({'Sun', 'Earth orbit', 'Jupiter orbit', 'Saturn orbit', ...
       'Earth position', 'Jupiter position', 'Saturn position'}, ...
       'FontSize', 14);
title('Jordan Mayer, HW 11 - Planet Positions at Jupiter Encounter');
mkdir('../img/allPos2');
save_all_views(fig_allPos2, '../img/allPos2/allPos2');

%%% Saturn encounter %%%
fprintf('\n---Saturn encounter---\n\n');

% get Saturn state vectors
[r_Saturn, v_Saturn] = state_vectors(mu_Sun, a_Saturn, e_Saturn, ...
    i_Saturn, raan_Saturn, aop_Saturn, t_0, M0_Saturn, t_3);
r_Saturn_3 = r_Saturn;  % for later
v_Saturn_3 = v_Saturn;  % for later

% plot Saturn position
fig_saturnPos3 = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Saturn_orbit;
plot_positions(r_Saturn);
legend({'Sun', 'Saturn orbit', 'Saturn position'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Saturn Position at Saturn Encounter');
mkdir('../img/saturnPos3');
save_all_views(fig_saturnPos3, '../img/saturnPos3/saturnPos3');

% plot all planet positions
[r_Earth, dummy] = state_vectors_quiet(mu_Sun, a_Earth, e_Earth, ...
    i_Earth, raan_Earth, aop_Earth, t_0, M0_Earth, t_3);
[r_Jupiter, dummy] = state_vectors_quiet(mu_Sun, a_Jupiter, e_Jupiter, ...
    i_Jupiter, raan_Jupiter, aop_Jupiter, t_0, M0_Jupiter, t_3);
fig_allPos3 = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_Jupiter_orbit;
plot_Saturn_orbit;
plot_positions([r_Earth, r_Jupiter, r_Saturn]);
legend({'Sun', 'Earth orbit', 'Jupiter orbit', 'Saturn orbit', ...
       'Earth position', 'Jupiter position', 'Saturn position'}, ...
       'FontSize', 14);
title('Jordan Mayer, HW 11 - Planet Positions at Saturn Encounter');
mkdir('../img/allPos3');
save_all_views(fig_allPos3, '../img/allPos3/allPos3');

%%%%% Part iii %%%%%

fprintf('\n-----Part iii-----\n');

%%% Earth to Jupiter %%%

fprintf('\n---Earth to Jupiter---\n\n');

r_1 = r_Earth_1  % km, x-y-z coords
r_1_mag = norm(r_1)  % km
r_2 = r_Jupiter_2  % km, x-y-z coords
r_2_mag = norm(r_2)  % km
TOF = t_2 - t_1  % days
TOF_s = TOF*24*60*60  % seconds

fprintf('\n');

rho1 = acos(dot(r_1,r_2)/(norm(r_1)*norm(r_2)))  % rad
rho2 = -rho1  % rad
rho3 = asin(norm(cross(r_1,r_2))/(norm(r_1)*norm(r_2)))  % rad
rho4 = pi - rho3  % rad
rho = rho4  % rad, chosen by consistency

% plot departure/arrival points to deduce type 1 or 2
fig_depArrE2J = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_Jupiter_orbit;
plot_positions([r_1 r_2]);
legend({'Sun', 'Earth orbit', 'Jupiter orbit', 'Earth launch', ...
    'Jupiter encounter'}, 'FontSize', 14);
ttl = ['Jordan Mayer, HW 11 - Departure/Arrival Points for ' ...
       'Earth-to-Jupiter Transfer'];
title(ttl);
mkdir('../img/depArrE2J');
save_all_views(fig_depArrE2J, '../img/depArrE2J/depArrE2J');

fprintf('\n');

c = norm(r_2 - r_1)  % km
s = 1/2 * (r_1_mag + r_2_mag + c)  % km

fprintf('\n');
TOF_par = 1/3 * sqrt(2/mu_Sun) * (s^(3/2) - (s - c)^(3/2))  % s
a_min = s/2  % km
alpha_min = pi  % rad
beta_min = 2*asin(sqrt((s-c)/(2*a_min)))  % rad
TOF_min = sqrt(a_min^3/mu_Sun) * ((alpha_min - sin(alpha_min)) - ...
                                  (beta_min - sin(beta_min)))  % s
                             
fprintf('\n');
                              
[a_e2j, alpha, beta] = get_sma_lambert_elliptic('1A', mu_Sun, TOF_s, c, ...
    s, [4.5e8 9.5e8])  % km, rad ,rad
p1 = 4*a_e2j*(s - r_1_mag)*(s - r_2_mag)/c^2 * (sin((alpha+beta)/2))^2  % km
p2 = 4*a_e2j*(s - r_1_mag)*(s - r_2_mag)/c^2 * (sin((alpha-beta)/2))^2  % km
p = p1  % km
e_e2j = sqrt(1 - p/a_e2j)

fprintf('\n');

h_hat = cross(r_1, r_2)/norm(cross(r_1, r_2))  % km^2/s
i_e2j = acos(h_hat(3))  % rad
i_deg = rad2deg(i_e2j)  % deg
raan1 = asin(h_hat(1)/sin(i_e2j))  % rad
raan2 = pi - raan1  % rad
raan3 = acos(-h_hat(2)/sin(i_e2j))  % rad
raan4 = -raan3  % rad
raan_e2j = raan3  % rad, chosen by consistency
raan_deg = rad2deg(raan_e2j)  % deg

fprintf('\n');

curlE = -mu_Sun/(2*a_e2j)  % km^2/s^2
IP = 2*pi*sqrt(a_e2j^3/mu_Sun)  % s
IP_day = IP/(60*24*24)  % days

fprintf('\n');

thetaStar_D1 = acos(1/e_e2j * (p/r_1_mag - 1))  % rad
thetaStar_D2 = -thetaStar_D1  % rad
thetaStar_A1 = acos(1/e_e2j * (p/r_2_mag - 1))  % rad
thetaStar_A2 = -thetaStar_A1  % rad
TA_11 = thetaStar_A1 - thetaStar_D1  % rad
TA_12 = thetaStar_A2 - thetaStar_D1  % rad
TA_21 = thetaStar_A1 - thetaStar_D2  % rad
TA_22 = thetaStar_A2 - thetaStar_D2  % rad
thetaStar_e2j_D = thetaStar_D1  % rad, chosen based on TA
thetaStar_D_deg = rad2deg(thetaStar_e2j_D)  % deg
thetaStar_e2j_A = thetaStar_A1  % rad, chosen based on TA
thetaStar_A_deg = rad2deg(thetaStar_e2j_A)  % deg

fprintf('\n');
r_D = r_1;  % km, x-y-z coords
r_A = r_2;  % km, x-y-z coords

fprintf('\n');

h = sqrt(mu_Sun * p)  % km^2/s
v_D_mag = sqrt(mu_Sun * (2/norm(r_D) - 1/a_e2j))  % km/s
v_A_mag = sqrt(mu_Sun * (2/norm(r_A) - 1/a_e2j))  % km/s
gamma_D = acos(h/(norm(r_D)*v_D_mag))  % rad
gamma_A = acos(h/(norm(r_A)*v_A_mag))  % rad
v_D_rth = [v_D_mag*sin(gamma_D); v_D_mag*cos(gamma_D); 0]
    % km/s, r-theta-h coords
v_A_rth = [v_A_mag*sin(gamma_A); v_A_mag*cos(gamma_A); 0]
    % km/s, r-theta-h coords

fprintf('\n');

r_hat_D = r_D/norm(r_D)  % x-y-z coords
theta_hat_D = cross(h_hat, r_hat_D)  % x-y-z coords
theta_D1 = asin(r_hat_D(3)/sin(i_e2j))  % rad
theta_D2 = pi - theta_D1  % rad
theta_D3 = acos(theta_hat_D(3)/sin(i_e2j))  % rad
theta_D4 = -theta_D3  % rad
theta_D = theta_D3  % rad, chosen for consistency
aop_e2j = theta_D - thetaStar_e2j_D  % rad
theta_A = aop_e2j + thetaStar_e2j_A  % rad

fprintf('\n');

v_D_e2j = rthetah_to_xyz(v_D_rth, i_e2j, raan_e2j, theta_D)
    % km/s, x-y-z coords
v_A_e2j = rthetah_to_xyz(v_A_rth, i_e2j, raan_e2j, theta_A)
    % km/s, x-y-z coords

%%% Jupiter to Saturn %%%

fprintf('\n---Jupiter to Saturn---\n\n');

r_D = r_Jupiter_2  % km, x-y-z coords
r_A = r_Saturn_3  % km, x-y-z coords
TOF = t_3 - t_2  % days
TOF_s = TOF*24*60*60  % seconds

fprintf('\n');

rho1 = acos(dot(r_D,r_A)/(norm(r_D)*norm(r_A)))  % rad
rho2 = -rho1
rho3 = asin(norm(cross(r_D,r_A))/(norm(r_D)*norm(r_A)))  % rad
rho4 = pi - rho3  % rad
rho = rho4  % rad

% plot departure/arrival points to deduce type 1 or 2
fig_depArrJ2S = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Jupiter_orbit;
plot_Saturn_orbit;
plot_positions([r_D r_A]);
legend({'Sun', 'Jupiter orbit', 'Saturn orbit', 'Jupiter encounter', ...
    'Saturn encounter'}, 'FontSize', 14);
ttl = ['Jordan Mayer, HW 11 - Departure/Arrival Points for ' ...
       'Jupiter-to-Saturn Transfer'];
title(ttl);
mkdir('../img/depArrJ2S');
save_all_views(fig_depArrJ2S, '../img/depArrJ2S/depArrJ2S');

fprintf('\n');

c = norm(r_A - r_D)  % km
s = 1/2 * (norm(r_D) + norm(r_A) + c)  % km
TOF_par = 1/3 * sqrt(2/mu_Sun) * (s^(3/2) - (s-c)^(3/2))  % s
a_min = s/2  % km
alpha_min = pi  % rad
beta_min = 2*asin(sqrt((s-c)/(2*a_min)))  % rad
TOF_min = sqrt(a_min^3 / mu_Sun) * (alpha_min - sin(alpha_min) - ...
                                    (beta_min - sin(beta_min)))  % s

fprintf('\n');

[a_j2s, alpha, beta] = get_sma_lambert_elliptic('1A', mu_Sun, TOF_s, c, ...
    s, [1.1e9 5e9])  % km, rad, rad
p1 = 4*a_j2s*(s - norm(r_D))*(s - norm(r_A))/c^2 * ...
    (sin((alpha+beta)/2))^2  % km
p2 = 4*a_j2s*(s - norm(r_D))*(s - norm(r_A))/c^2 * ...
    (sin((alpha-beta)/2))^2  % km
p = p1  % km, chosen by type
e_j2s = sqrt(1 - p/a_j2s)

fprintf('\n');

h_hat = cross(r_D,r_A)/norm(cross(r_D,r_A))  % x-y-z coords
i_j2s = acos(h_hat(3))  % rad
i_j2s_deg = rad2deg(i_j2s)  % deg
raan1 = asin(h_hat(1)/sin(i_j2s))  % rad
raan2 = pi - raan1  % rad
raan3 = acos(-h_hat(2)/sin(i_j2s))  % rad
raan4 = -raan3  % rad
raan_j2s = (raan4)  % rad, between 0 and 2pi
raan_j2s_deg = rad2deg(raan_j2s)  % deg

fprintf('\n');

curlE = -mu_Sun/(2*a_j2s)  % km^2/s^2
IP = 2*pi*sqrt(a_j2s^3/mu_Sun)  % s
IP_yr = IP/24/60/60/365  % years

fprintf('\n');

thetaStar_D1 = acos(1/e_j2s * (p/norm(r_D) - 1))  % rad
thetaStar_D2 = -thetaStar_D1  % rad
thetaStar_A1 = acos(1/e_j2s * (p/norm(r_A) - 1))  % rad
thetaStar_A2 = -thetaStar_A1  % rad
TA_11 = thetaStar_A1 - thetaStar_D1  % rad
TA_12 = thetaStar_A2 - thetaStar_D1  % rad
TA_21 = thetaStar_A1 - thetaStar_D2  % rad
TA_22 = thetaStar_A2 - thetaStar_D2  % rad
thetaStar_j2s_D = thetaStar_D2  % rad
thetaStar_j2s_D_deg = rad2deg(thetaStar_j2s_D)  % deg
thetaStar_j2s_A = thetaStar_A1  % rad
thetaStar_j2s_A_deg = rad2deg(thetaStar_j2s_A)  % deg

fprintf('\n');

h = sqrt(mu_Sun*p)  % km^2/s
v_D_mag = sqrt(mu_Sun*(2/norm(r_D) - 1/a_j2s))  % km/s
v_A_mag = sqrt(mu_Sun*(2/norm(r_A) - 1/a_j2s))  % km/s
gamma_D = -acos(h/(norm(r_D) * v_D_mag))  % rad
gamma_A = acos(h/(norm(r_A) * v_A_mag))  % rad
v_D_rth = [v_D_mag*sin(gamma_D); v_D_mag*cos(gamma_D); 0]
  % km/s, r-theta-h coords
v_A_rth = [v_A_mag*sin(gamma_A); v_A_mag*cos(gamma_A); 0]
  % km/s, r-theta-h coords
  
fprintf('\n');

r_hat_D = r_D/norm(r_D)  % x-y-z coords
theta_hat_D = cross(h_hat, r_hat_D)  % x-y-z coords
theta_D1 = asin(r_hat_D(3)/sin(i_j2s))  % rad
theta_D2 = pi - theta_D1  % rad
theta_D3 = acos(theta_hat_D(3)/sin(i_j2s))  % rad
theta_D4 = -theta_D3  % rad
theta_D = theta_D4  % rad
aop_j2s = theta_D - thetaStar_j2s_D  % rad
theta_A = thetaStar_j2s_A + aop_j2s  % rad

fprintf('\n');

v_D_j2s = rthetah_to_xyz(v_D_rth, i_j2s, raan_j2s, theta_D)
    % km, x-y-z coords
v_A_j2s = rthetah_to_xyz(v_A_rth, i_j2s, raan_j2s, theta_A)
    % km, x-y-z coords

%%%%% Part iv %%%%%

axes_limits = [-norm(r_Saturn_3)*1.1, norm(r_Saturn_3)*1.1];

% plot Earth to Jupiter transfer
% plot planet orbits
fig_e2j = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_Jupiter_orbit;
% plot departure/arrival points
plot_positions([r_Earth_1, r_Jupiter_2])
% plot transfer arc
plot_orbit_3D(a_e2j, e_e2j, i_e2j, raan_e2j, aop_e2j, ...
    [thetaStar_e2j_D, thetaStar_e2j_A]);
xlim(axes_limits);
ylim(axes_limits);
zlim(axes_limits);
legend({'Sun', 'Earth orbit', 'Jupiter orbit', ...
    'Earth launch (departure)', 'Jupiter encounter (arrival)', ...
    'Transfer arc'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Earth to Jupiter Transfer');
mkdir('../img/e2j');
save_all_views(fig_e2j, '../img/e2j/e2j');

% plot Jupiter to Saturn transfer
% plot planet orbits
fig_j2s = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Jupiter_orbit;
plot_Saturn_orbit;
% plot departure/arrival points
plot_positions([r_Jupiter_2, r_Saturn_3]);
% plot transfer arc
plot_orbit_3D(a_j2s, e_j2s, i_j2s, raan_j2s, aop_j2s, ...
    [thetaStar_j2s_D, thetaStar_j2s_A]);
xlim(axes_limits);
ylim(axes_limits);
zlim(axes_limits);
legend({'Sun', 'Jupiter orbit', 'Saturn orbit', ...
    'Jupiter encounter (departure)', 'Saturn encounter (arrival)', ...
    'Transfer arc'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Jupiter to Saturn Transfer');
mkdir('../img/j2s');
save_all_views(fig_j2s, '../img/j2s/j2s');

% plot both transfers
% plot planet orbits
fig_e2j2s = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Earth_orbit;
plot_Jupiter_orbit;
plot_Saturn_orbit;
% plot departure/arrival points
plot_positions([r_Earth_1, r_Jupiter_2, r_Saturn_3]);
% plot transfer arcs
plot_orbit_3D(a_e2j, e_e2j, i_e2j, raan_e2j, aop_e2j, ...
    [thetaStar_e2j_D, thetaStar_e2j_A]);
plot_orbit_3D(a_j2s, e_j2s, i_j2s, raan_j2s, aop_j2s, ...
    [thetaStar_j2s_D, thetaStar_j2s_A]);
xlim(axes_limits);
ylim(axes_limits);
zlim(axes_limits);
legend({'Sun', 'Earth orbit', 'Jupiter orbit', 'Saturn orbit', ...
    'Earth launch', 'Jupiter encounter', 'Saturn encounter', ...
    'Earth to Jupiter transfer arc', ...
    'Jupiter to Saturn transfer arc'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Both Transfers');
mkdir('../img/e2j2s');
save_all_views(fig_e2j2s, '../img/e2j2s/e2j2s');

%%%%%%%% Part b: Launch %%%%%%%%

fprintf('\n-------Part b: Launch--------\n\n');

v_infEarth_plus = v_D_e2j - v_Earth_1  % km/s, x-y-z coords
v_infEarth_mag = norm(v_infEarth_plus)  % km/s

fprintf('\n');

r_0 = R_Earth + 300  % km
v_p = sqrt(2 * (v_infEarth_mag^2 / 2 + mu_Earth/r_0))  % km/s
v_0 = sqrt(mu_Earth/r_0)  % km/s
deltaV = v_p - v_0  % km/s

%%%%%%%% Part c: Jupiter encounter %%%%%%%%

fprintf('\n--------Part c: Jupiter encounter--------\n');

%%%%% Part i %%%%%

fprintf('\n-----Part i-----\n\n');

v_A1 = v_A_e2j  % km/s, x-y-z coords
v_D2 = v_D_j2s  % km/s, x-y-z coords
deltaV_eq = v_D2 - v_A1  % km/s, x-y-z coords
deltaV_eq_mag = norm(deltaV_eq)  % km/s

%%%%% Part ii %%%%%

fprintf('\n-----Part ii-----\n\n');

v_infJupiter_minus = v_A1 - v_Jupiter_2  % km/s, x-y-z coords
v_infJupiter_minus_mag = norm(v_infJupiter_minus)  % km/s
v_infJupiter_plus = v_D2 - v_Jupiter_2  % km/s, x-y-z coords
v_infJupiter_plus_mag = norm(v_infJupiter_plus)  % km/s
v_infJupiter_diff = abs(v_infJupiter_minus_mag - v_infJupiter_plus_mag)
    % km/s

fprintf('\n');

v_infJupiter = 1/2 * (v_infJupiter_minus_mag + v_infJupiter_plus_mag)
  % km/s
abs_a = mu_Jupiter/v_infJupiter^2  % km
del1 = asin(norm(cross(v_infJupiter_minus, v_infJupiter_plus)) / ...
           (v_infJupiter_minus_mag * v_infJupiter_plus_mag))  % rad
del2 = pi - del1  % rad
del3 = acos(dot(v_infJupiter_minus, v_infJupiter_plus) / ...
            (v_infJupiter_minus_mag * v_infJupiter_plus_mag))  % rad
del4 = -del3  % rad
del = del3  % rad, chosen by consistency

fprintf('\n');

e = 1/sin(del/2)
r_p = abs_a * (e - 1)  % km
h_p = r_p - R_Jupiter  % km
h_p_act = 42828;  % km
h_p_diff = abs(h_p - h_p_act)  % km

fprintf('\n');

psi1 = asin(norm(cross(v_A1, v_infJupiter_minus)) / ...
            (norm(v_A1) * v_infJupiter_minus_mag))  % rad
psi2 = asin(norm(cross(v_D2, v_A1)) / (norm(v_D2) * norm(v_A1)))  % rad
psi3 = asin(norm(cross(v_infJupiter_plus, v_D2)) / ...
            (v_infJupiter_plus_mag * norm(v_D2)))  % rad
psi_sum = psi1 + psi2 + psi3  % rad

%%%%%%%% Part d: Saturn Encounter %%%%%%%%

fprintf('\n--------Part d: Saturn Encounter--------\n');

%%%%% Part i %%%%%

fprintf('\n-----Part i-----\n\n');

h_p_act = 20900;  % km
r_p = h_p_act + R_Saturn  % km

fprintf('\n');

v_Saturn = v_Saturn_3  % km/s, x-y-z coords
v_minus = v_A_j2s  % km/s, x-y-z coords
v_infSaturn_minus = v_minus - v_Saturn  % km/s, x-y-z coords
v_infSaturn = norm(v_infSaturn_minus)  % km/s

fprintf('\n');

abs_a = mu_Saturn/v_infSaturn^2  % km
e = r_p/abs_a + 1
del = 2*asin(1/e)  % rad

fprintf('\n');

v_infSaturn_plus = Saturn_flyby(r_p, v_infSaturn_minus, v_Saturn, del)
    % km/s, x-y-z frame, both options
v_plus = [norm(v_infSaturn_plus(:,1) + v_Saturn), ...
          norm(v_infSaturn_plus(:,2) + v_Saturn)]  % km/s, both options
v_infSaturn_plus = v_infSaturn_plus(:,2)  % km/s, x-y-z coords
v_plus = v_infSaturn_plus + v_Saturn  % km/s, x-y-z coords
v_plus_mag = norm(v_plus)  % km/s

fprintf('\n');

r_D = r_Saturn_3  % km, x-y-z coords
v_D = v_plus  % km/s, x-y-z coords
curlE = norm(v_D)^2 / 2 - mu_Sun/norm(r_D)  % km^2/s^2
a_s2n = -mu_Sun/(2*curlE)  % km
h = norm(cross(r_D, v_D))  % km^2/s
p = h^2 / mu_Sun  % km
e_s2n = sqrt(1 - p/a_s2n)

fprintf('\n');

h_hat = cross(r_D, v_D)/h  % x-y-z coords
i_s2n = acos(h_hat(3))  % rad
i_s2n_deg = rad2deg(i_s2n)  % deg
raan_s2n1 = asin(h_hat(1)/sin(i_s2n))  % rad
raan_s2n2 = pi - raan_s2n1  % rad
raan_s2n3 = acos(-h_hat(2)/sin(i_s2n))  % rad
raan_s2n4 = -raan_s2n3  % rad
raan_s2n = raan_s2n3  % rad, chosen for consistency
raan_s2n_deg = rad2deg(raan_s2n)  % deg

fprintf('\n');

r_hat_D = r_D/norm(r_D)  % x-y-z coords
theta_hat_D = cross(h_hat, r_hat_D)  % x-y-z coords
theta_D1 = asin(r_hat_D(3)/sin(i_s2n))  % rad
theta_D2 = pi - theta_D1  % rad
theta_D3 = acos(theta_hat_D(3)/sin(i_s2n))  % rad
theta_D4 = -theta_D3  % rad
theta_D = theta_D1  % rad, chosen by consistency
thetaStar_D1 = acos(1/e_s2n * (p/norm(r_D) - 1))  % rad
thetaStar_D2 = -thetaStar_D1  % rad
v_D_r = dot(v_D, r_D)/norm(r_D)  % km/s
thetaStar_s2n_D = thetaStar_D2  % rad, chosen b/c descending
thetaStar_D_deg = rad2deg(thetaStar_s2n_D)  % deg
aop_s2n = theta_D - thetaStar_s2n_D  % rad

% plot final trajectory
axes_limits = [-a_Neptune * 1.2, a_Neptune * 1.2];
thetaStar_s2n_inf = acos(-1/e);
thetaStar_s2n_stop = thetaStar_s2n_inf - 0.8;
fig_s2n = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
plot_Saturn_orbit;
plot_Neptune_orbit;
plot_positions(r_D);
plot_orbit_3D(a_s2n, e_s2n, i_s2n, raan_s2n, aop_s2n, ...
    [thetaStar_s2n_D, thetaStar_s2n_stop]);
xlim(axes_limits);
ylim(axes_limits);
zlim(axes_limits);
legend({'Sun', 'Saturn orbit', 'Neptune orbit', 'Saturn encounter', ...
    'Final trajectory'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Final Trajectory following Saturn Encounter');
mkdir('../img/s2n');
save_all_views(fig_s2n, '../img/s2n/s2n');

%%%%% Part ii %%%%%

fprintf('\n-----Part ii-----\n\n');

thetaStar_C = acos(1/e_s2n * (p/a_Neptune - 1))  % rad
H_C = 2*atanh(((e_s2n - 1)/(e_s2n + 1))^(1/2) * tan(thetaStar_C/2))  % rad
H_D = 2*atanh(((e_s2n - 1)/(e_s2n + 1))^(1/2) * tan(thetaStar_s2n_D/2))
    % rad
N_C = e_s2n*sinh(H_C) - H_C
N_D = e_s2n*sinh(H_D) - H_D
deltaN = N_C - N_D  % rad
deltaT = sqrt(abs(a_s2n)^3/mu_Sun) * deltaN  % s
deltaT_day = deltaT/24/60/60  % days
t_C = t_3 + deltaT_day  % JD

fprintf('\n');

r_C = rthetah_to_xyz([a_Neptune 0 0]', i_s2n, raan_s2n, ...
    aop_s2n + thetaStar_C)  % km, x-y-z coords
r_Earth_C = state_vectors_quiet(mu_Sun, a_Earth, e_Earth, i_Earth, ...
    raan_Earth, aop_Earth, t_0, M0_Earth, t_C)  % km, x-y-z coords
dist_Earth_C = norm(r_C - r_Earth_C)  % km
dist_Earth_C_AU = dist_Earth_C/149597900  % AU

%%%%%%%% Part e %%%%%%%%
axes_limits = [-a_Neptune * 1.2, a_Neptune * 1.2];
fig_completePath = figure('Position', fig_position);
plot_positions([0 0 0]');
hold on;
% plot planet orbits
plot_Earth_orbit;
plot_Jupiter_orbit;
plot_Saturn_orbit;
plot_Neptune_orbit;
% plot all departure/arrival points
plot_positions([r_Earth_1, r_Jupiter_2, r_Saturn_3]);
% plot all trajectories
plot_orbit_3D(a_e2j, e_e2j, i_e2j, raan_e2j, aop_e2j, ...
    [thetaStar_e2j_D, thetaStar_e2j_A]);
plot_orbit_3D(a_j2s, e_j2s, i_j2s, raan_j2s, aop_j2s, ...
    [thetaStar_j2s_D, thetaStar_j2s_A]);
plot_orbit_3D(a_s2n, e_s2n, i_s2n, raan_s2n, aop_s2n, ...
    [thetaStar_s2n_D, thetaStar_s2n_stop]);
xlim(axes_limits);
ylim(axes_limits);
zlim(axes_limits);
legend({'Sun', 'Earth orbit', 'Jupiter orbit', 'Saturn orbit', ...
    'Neptune orbit', 'Launch', 'Jupiter encounter', 'Saturn encounter', ...
    'Earth to Jupiter trajectory', 'Jupiter to Saturn trajectory', ...
    'Post-Saturn trajectory'}, 'FontSize', 14);
title('Jordan Mayer, HW 11 - Complete Heliocentric Path');
mkdir('../img/completePath');
save_all_views(fig_completePath, '../img/completePath/completePath');