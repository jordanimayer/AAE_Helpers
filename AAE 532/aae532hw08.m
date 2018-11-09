%%%%%
% Jordan Mayer
% AAE 532
% HW 08
%
% Perform various arithmetic and create plots to perform orbital analysis.
%%%%%

%%% Preparations %%%

clear all; close all; format compact; format long;

% load necessary constants (see load_constants for explanations and units)
global mu_Sun ...
       R_Moon mu_Moon r_Moon ...
       R_Earth mu_Earth r_Earth ...
       R_Jupiter mu_Jupiter r_Jupiter ...
       R_Pluto mu_Pluto r_Pluto
load_constants;

%%% Problem 1 %%%

fprintf('\nProblem 1:\n');

% Part a

fprintf('\nPart a:\n\n');

V_pT = 41.5874;  % km/s
V_Earth = 29.7794;  % km/s
V_infEarth = V_pT - V_Earth  % km/s
r_i = R_Earth + 180  % km
V_1_plus = sqrt(2 * (V_infEarth^2 / 2 + mu_Earth/r_i))  % km/s
V_1_minus = sqrt(mu_Earth/r_i)  % km/s
deltaV_1 = V_1_plus - V_1_minus  % km/s

fprintf('\n');

V_aT = 1.0605;  % km/s
V_Pluto = 4.7556;  % km/s
V_infPluto = V_Pluto - V_aT  % km/s
r_f = 7 * R_Pluto  % km
V_2_minus = sqrt(2 * (V_infPluto^2 / 2 + mu_Pluto/r_f))  % km/s
V_2_plus = sqrt(mu_Pluto/r_f)  % km/s
deltaV_2 = V_2_minus - V_2_plus  % km/s

fprintf('\n');

deltaV_Tot = deltaV_1 + deltaV_2  % km/s

% Part b

fprintf('\nPart b:\n\n');

r_af = 30*R_Pluto  % km
V_2_minus = sqrt(2 * (V_infPluto^2 / 2 + mu_Pluto/r_af))  % km/s
a_f = 16*R_Pluto  % km
V_2_plus = sqrt(2 * mu_Pluto * (1/r_af - 1/(2*a_f)))  % km/s
deltaV_2 = V_2_minus - V_2_plus  % km/s

%%% Problem 2 %%%

fprintf('\nProblem 2:\n');

% Part a

fprintf('\nPart a:\n\n');

t_launch = 2453754.5;  % Julian Days, UT
t_e = 2454159.5;  % Julian Days, UT
deltaT = t_e - t_launch  % Julian Days

% Part b

fprintf('\nPart b:\n\n');

e_T = 0.9503;
a_T = 3.00888e9;  % km
p_T = 2.91857e8;  % km

r_e = r_Jupiter  % km
thetaStar_e_1 = acosd((p_T/r_e - 1)/e_T)  % deg, option 1
thetaStar_e_2 = -acosd((p_T/r_e - 1)/e_T)  % deg, option 2
thetaStar_e = thetaStar_e_1  % deg, chosen b/c ascending
n_T = sqrt(mu_Sun/a_T^3)  % rad/s
E_e_1 = 2*atand(((1 - e_T)/(1 + e_T))^(1/2) * tand(thetaStar_e / 2))
  % deg, option 1
E_e_2 = 180 + E_e_1  % deg, option 2
E_e = E_e_1  % deg, chosen b/c ascending
E_e_rad = deg2rad(E_e)  % rad
M_e = E_e_rad - e_T * sind(E_e)
TOF = M_e/n_T  % s
TOF_day = TOF/60/60/24  % days
phi = deg2rad(thetaStar_e) - n_T * TOF  % rad
phi_deg = rad2deg(phi)

% Part d

fprintf('\nPart d:\n\n');

deltaV_1 = 8.35900;  % km/s
r_i = R_Earth + 180  % km
V_esc = sqrt(2*mu_Earth/r_i)  % km/s
V_1_minus = 7.79612;  % km/s
deltaV_esc = V_esc - V_1_minus  % km/s
diff_deltaV_1_esc = deltaV_1 - deltaV_esc  % km/s

% Part e

fprintf('\nPart e:\n\n');

r_minus = r_Jupiter  % km
curlE_T = -mu_Sun/(2*a_T)  % km^2/s^2
V_minus = sqrt(2 * (curlE_T + mu_Sun/r_minus))  % km/s
gamma_minus_1 = acosd(sqrt(mu_Sun * p_T)/(r_minus * V_minus))
  % deg, option 1
gamma_minus_2 = -gamma_minus_1  % deg, option 2
gamma_minus = gamma_minus_1  % deg, chosen b/c ascending
thetaStar_minus = 131.153;  % deg
V_minus_r = V_minus * sind(gamma_minus)  % km/s
V_minus_theta = V_minus * cosd(gamma_minus)  % km/s

% Part f

fprintf('\nPart f:\n\n');

r_plus = r_minus;  % km
V_Jupiter = sqrt(mu_Sun/r_Jupiter)  % km/s
V_infJupiter = sqrt(V_minus^2 + V_Jupiter^2 - ...
                    2*V_minus*V_Jupiter*cosd(gamma_minus))  % km/s

fprintf('\n');

xi_1 = asind(V_minus/V_infJupiter * sind(gamma_minus))  % deg, option 1
xi_2 = 180 - xi_1  % deg, option 2
psi_1 = asind(V_Jupiter/V_infJupiter * sind(gamma_minus))  % deg, option 1
psi_2 = 180 - psi_1  % deg, option 2
sum_11 = gamma_minus + xi_1 + psi_1
sum_12 = gamma_minus + xi_1 + psi_2
sum_21 = gamma_minus + xi_2 + psi_1
sum_22 = gamma_minus + xi_2 + psi_2
xi = xi_1  % deg, chosen b/c sum = 180
psi = psi_1  % deg, chosen b/c sum = 180

fprintf('\n');

abs_a = mu_Jupiter/V_infJupiter^2  % km
r_p = 32*R_Jupiter  % km
e = r_p/abs_a + 1
del = 2*asind(1/e)  % deg
nu = (180 - del)/2  % deg
deltaV_eq = sind(del)/sind(nu) * V_infJupiter  % km/s
alpha = 180 - psi - nu  % deg

fprintf('\n');

eta = xi + del  % deg
V_plus = sqrt(V_Jupiter^2 + V_infJupiter^2 - ...
              2*V_Jupiter*V_infJupiter*cosd(eta))  % km/s
gamma_plus_1 = asind(V_infJupiter/V_plus * sind(eta))  % deg, option 1
gamma_plus_2 = 180 - gamma_plus_1  % deg, option 2
gamma_plus = gamma_plus_1  % deg, chosen b/c ascending

fprintf('\n');

curlE_plus = V_plus^2 / 2 - mu_Sun/r_plus  % km^2/s^2
a_plus = -mu_Sun/(2*curlE_plus)  % km
h_plus = r_plus*V_plus*cosd(gamma_plus)  % km^2/s
p_plus = h_plus^2/mu_Sun  % km
e_plus = sqrt(1 - p_plus/a_plus)
thetaStar_plus_1 = acosd(1/e_plus * (p_plus/r_plus - 1))  % deg, option 1
thetaStar_plus_2 = -thetaStar_plus_1  % deg, option 2
thetaStar_plus = thetaStar_plus_1  % deg, chosen b/c ascending
H_plus = 2*atanh(((e_plus - 1)/(e_plus + 1))^(1/2) * ...
                 tand(thetaStar_plus/2))
N_plus = e_plus * sinh(H_plus) - H_plus
delta_aop = thetaStar_minus - thetaStar_plus  % deg

fprintf('\n');

V_plus_r = V_plus*sind(gamma_plus)  % km/s
V_plus_theta = V_plus*cosd(gamma_plus)  % km/s

% Part g

fprintf('\nPart g:\n\n');

thetaStar_inf_plus = acosd(-1/e_plus)  % deg

figure(); hold on;

% plot old orbit
plotOrbit2D(e_T, p_T, [0 360], 0, '');

% plot new orbit
ttl = 'Jordan Mayer - HW 08, Problem 2g, Heliocentric View';
plotOrbit2D(e_plus, p_plus, [-140 140], delta_aop, ttl);

legend('Before Jupiter encounter', 'After Jupiter encounter');
axis(10^9 * [-7 2 -4.5 4.5]);
axis square;

%%% Problem 3 %%%

fprintf('\nProblem 3:\n');

% Part 1

fprintf('\nPart 1:\n');

% Part a

fprintf('\nPart a:\n\n');

r_i = R_Earth + 200  % km
V_1_minus = sqrt(mu_Earth/r_i)  % km/s
r_p_T = r_i  % km
r_a_T = r_Moon  % km
a_T = (r_p_T + r_a_T)/2  % km
curlE_T = -mu_Earth/(2*a_T)  % km^2/s^2
V_1_plus = sqrt(2 * (curlE_T + mu_Earth/r_i))  % km/s
deltaV_1 = V_1_plus - V_1_minus  % km/s

fprintf('\n');

V_minus = sqrt(2 * (curlE_T + mu_Earth/r_Moon))  % km/s
V_Moon = sqrt(mu_Earth/r_Moon)  % km/s
V_infMoon = V_Moon - V_minus  % km/s

fprintf('\n');

r_f = R_Moon + 140  % km
V_2_plus = sqrt(mu_Moon/r_f)  % km/s
curlE = V_infMoon^2 / 2  % km^2/s^2
V_2_minus = sqrt(2 * (curlE + mu_Moon/r_f))  % km/s
deltaV_2 = V_2_minus - V_2_plus  % km/s

fprintf('\n');

deltaV_Tot = deltaV_1 + deltaV_2  % km/s

fprintf('\n');

IP_T = 2*pi*sqrt(a_T^3 / mu_Earth)  % s
IP_T_day = IP_T/60/60/24  % days
TOF = IP_T/2  % s
TOF_day = IP_T_day/2  % days

fprintf('\n');

n_T = sqrt(mu_Earth / r_Moon^3)  % rad/s
phi = pi - n_T*TOF  % rad
phi_deg = rad2deg(phi)  % deg

fprintf('\n');

h_T = r_Moon * V_minus  % km^2/s
p_T = h_T^2/mu_Earth  % km
e_T = sqrt(1 - p_T/a_T)

% Part b

fprintf('\nPart b:\n\n');

r_p = r_f;  % km
abs_a = mu_Moon/V_infMoon^2  % km
e = r_p/abs_a + 1
del = 2*asind(1/e)  % deg
V_plus = sqrt(V_infMoon^2 + V_Moon^2 - 2*V_infMoon*V_Moon*cosd(del))
  % km/s
r_plus = r_Moon;  % km
gamma_plus_1 = asind(V_infMoon/V_plus * sind(del))  % deg, option 1
gamma_plus_2 = 180 - gamma_plus_1  % deg, option 2
gamma_plus = gamma_plus_1  % deg, chosen b/c must be b/w 0 and 90

fprintf('\n');

curlE_plus = V_plus^2 / 2 - mu_Earth/r_plus  % km^2/s^2
a_plus = -mu_Earth/(2*curlE_plus)  % km
h_plus = r_plus*V_plus*cosd(gamma_plus)  % km^2/s
p_plus = h_plus^2 / mu_Earth  % km
e_plus = sqrt(1 - p_plus/a_plus)
thetaStar_plus = acosd(1/e_plus * (p_plus/r_plus - 1))  % deg
delta_aop = 180 - thetaStar_plus  % deg
r_p_plus = a_plus * (1 - e_plus)  % km

% Part c

fprintf('\nPart c:\n\n');

deltaV_eq = sqrt(2*V_infMoon^2 * (1 - cosd(del)))  % km/s
alpha = (180 - del)/2  % deg
deltaV_eq_v = deltaV_eq * cosd(alpha)  % km/s
deltaV_eq_c = deltaV_eq * sind(alpha)  % km/s

% Part 2

fprintf('\nPart 2:\n');

% Part d

fprintf('\nPart d:\n\n');

r_minus = r_Moon  % km
r_p_minus = r_i  % km
thetaStar_minus = 172.5;  % deg
e_minus = (r_minus - r_p_minus)/(r_p_minus - r_minus*cosd(thetaStar_minus))
a_minus = r_p_minus/(1 - e_minus)  % km
r_a_minus = a_minus * (1 + e_minus)  % km
IP_minus = 2*pi*sqrt(a_minus^3 / mu_Earth)  % s
IP_minus_day = IP_minus/60/60/24  % days
curlE_minus = -mu_Earth/(2*a_minus)  % km^2/s^2
V_minus = sqrt(2 * (curlE_minus + mu_Earth/r_minus))  % km/s
p_minus = r_p_minus * (1 - e_minus)  % km
h_minus = sqrt(mu_Earth * p_minus)  % km^2/s
gamma_minus = acosd(h_minus/(r_minus * V_minus))  % deg
V_infMoon = sqrt(V_Moon^2 + V_minus^2 - 2*V_Moon*V_minus*cosd(gamma_minus))
  % km/s
E_minus = 2*atand(((1 - e_minus)/(1 + e_minus))^(1/2) * ...
                  tand(thetaStar_minus/2))  % deg
E_minus_rad = deg2rad(E_minus)  % rad
n_minus = sqrt(mu_Earth/a_minus^3)  % rad/s
M_minus = E_minus_rad - e_minus * sind(E_minus)
TOF_minus = M_minus/n_minus  % s
TOF_minus_day = TOF_minus/60/60/24  % days

% Part e

fprintf('\nPart e:\n\n')

eta = asind(V_minus/V_infMoon * sind(gamma_minus))  % deg
del = 2*eta  % deg
deltaV_eq = sqrt(2 * V_infMoon^2 * (1 - cosd(del)))  % km/s
beta = -acosd(deltaV_eq/(2 * V_minus))  % deg
alpha = 180 - beta  % deg
alpha = alpha - 360  % deg
e = 1/sind(eta)
abs_a = mu_Moon/V_infMoon^2  % km
r_p = abs_a * (e - 1)  % km
h_p = r_p - R_Moon  % km, altitude

fprintf('\n');

V_plus = V_minus  % km/s
gamma_plus = -gamma_minus  % deg
r_plus = r_minus  % km
r_p_plus = r_p_minus  % km
curlE_plus = (V_plus^2 / 2 - mu_Earth/r_plus)  % km^2/s^2
a_plus = -mu_Earth/(2*curlE_plus)  % km
h_plus = r_plus*V_plus*cosd(gamma_plus)  % km^2/s
p_plus = h_plus^2/mu_Earth  % km
e_plus = sqrt(1 - p_plus/a_plus)
thetaStar_plus = -acosd(1/e_plus * (p_plus/r_plus - 1))  % deg
IP_plus = 2*pi*sqrt(a_plus^3 / mu_Earth)  % s
IP_plus_day = IP_plus/60/60/24  % days
r_a_plus = a_plus * (1 + e_plus)  % km
delta_aop = thetaStar_minus - thetaStar_plus - 360  % deg

% Part g

fprintf('\nPart g:\n\n');

deltaV_eq_v = -deltaV_eq * cosd(abs(beta))  % km/s
deltaV_eq_c = -deltaV_eq * sind(abs(beta))  % km/s