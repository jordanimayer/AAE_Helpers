%%%%%
% Jordan Mayer
% AAE 532
% HW 09
%
% Perform various arithmetic to analyze orbits.
%%%%%

%%% Preliminary setup %%%

clear all; close all; format compact; format long;
global mu_Sun ...
       r_Mercury ...
       r_Venus R_Venus mu_Venus ...
       r_Earth ...
       R_Moon mu_Moon ...
       R_Saturn mu_Saturn;
   % see load_constants for descriptions and units
load_constants;

%%% Problem 1 %%%

fprintf('\nProblem 1:\n');

% Part a

fprintf('\nPart a:\n');

% Option i

fprintf('\nOption i:\n\n');

r_c = R_Moon + 100  % km
V_c = sqrt(mu_Moon/r_c)  % km/s
deltaV = 2*V_c*sind(45)  % km/s

% Option ii

fprintf('\nOption ii:\n\n');

V_1_minus = V_c  % km/s
r_a = 25*R_Moon  % km
r_p = r_c  % km
a = 1/2 * (r_a + r_p)  % km
curlE = -mu_Moon/(2*a)  % km^2/s^2
V_1_plus = sqrt(2 * (curlE + mu_Moon/r_p))  % km/s
deltaV_1 = V_1_plus - V_1_minus  % km/s

fprintf('\n')

V_a = sqrt(2 * (curlE + mu_Moon/r_a))  % km/s
deltaV_2 = 2*V_a*sind(45)  % km/s

fprintf('\n')

V_3_minus = V_1_plus  % km/s
V_3_plus = V_c  % km/s
deltaV_3 = V_3_minus - V_3_plus  % km/s

fprintf('\n')

deltaV = deltaV_1 + deltaV_2 + deltaV_3  % km/s

fprintf('\n')

TOF = 2*pi*sqrt(a^3/mu_Moon)  % s
TOF_day = TOF/60/60/24  % days

% Part b

fprintf('\nPart b:\n\n');

V_1_minus = V_c  % km/s
r_a = 17000  % km
r_p = r_c  % km
a = 1/2 * (r_a + r_p)  % km
curlE = -mu_Moon/(2*a)  % km^2/s^2
V_1_plus = sqrt(2 * (curlE + mu_Moon/r_p))  % km/s
deltaV_1 = V_1_plus - V_1_minus  % km/s

fprintf('\n')

V_a = sqrt(2 * (curlE + mu_Moon/r_a))  % km/s
deltaV_2 = 2*V_a*sind(45)  % km/s

fprintf('\n')

V_3_minus = V_1_plus  % km/s
V_3_plus = V_c  % km/s
deltaV_3 = V_3_minus - V_3_plus  % km/s

fprintf('\n')

deltaV = deltaV_1 + deltaV_2 + deltaV_3  % km/s

fprintf('\n')

TOF = 2*pi*sqrt(a^3/mu_Moon)  % s
TOF_day = TOF/60/60/24  % days

%%% Problem 2 %%%

fprintf('\nPart 2:\n');

% Part a

fprintf('\nPart a:\n\n');

r_p_minus = 8*R_Saturn  % km
r_a_minus = 32*R_Saturn  % km
a_Titan = 20*R_Saturn  % km
e_Titan = 0.2;
mu_Titan = 8825;  % km^3/s^2
R_Titan = 2575;  % km

fprintf('\n')

p_Titan = a_Titan * (1 - e_Titan^2)  % km
e_minus = 1 - 2*r_p_minus/(r_p_minus + r_a_minus)
p_minus = r_p_minus * (1 + e_minus)  % km
thetaStar_minus = acosd((p_Titan/p_minus - 1) * ...
                        (e_Titan - p_Titan/p_minus * e_minus)^(-1))
  % deg
r = p_minus/(1 + e_minus * cosd(thetaStar_minus))  % km

fprintf('\n');

curlE_Titan = -mu_Saturn/(2*a_Titan)  % km^2/s^2
V_Titan = sqrt(2*(curlE_Titan + mu_Saturn/r))  % km/s
a_minus = 1/2 * (r_p_minus + r_a_minus)  % km
curlE_minus = -mu_Saturn/(2*a_minus)  % km^2/s^2
V_minus = sqrt(2*(curlE_minus + mu_Saturn/r))  % km/s

fprintf('\n');

h_Titan = sqrt(mu_Saturn*p_Titan)  % km^2/s
gamma_Titan = acosd(h_Titan/(r*V_Titan))  % deg
h_minus = sqrt(mu_Saturn*p_minus)  % km^2/s
gamma_minus = acosd(h_minus/(r*V_minus))  % deg

fprintf('\n');

xi_minus = gamma_minus - gamma_Titan  % deg
V_infTitan = sqrt(2 * V_minus^2 * (1 - cosd(xi_minus)))  % km/s
eta_minus = 1/2 * (180 - xi_minus)  % deg
psi_minus = eta_minus  % deg

fprintf('\n');

curlE = V_infTitan^2 / 2  % km^2/s^2
abs_a = mu_Titan/(2*curlE)  % km
r_p = R_Titan + 1000  % km
e = r_p/abs_a + 1
del = 2*asind(1/e)  % deg

fprintf('\n');

zeta = 1/2 * (180 - del)  % deg
deltaV_eq = sind(del)/sind(zeta) * V_infTitan  % km/s

fprintf('\n');

beta = zeta + psi_minus  % deg
alpha = 180 - beta  % deg

fprintf('\n');

eta_plus = 360 - del - eta_minus  % deg, > 180 indicates bad diagram
eta_plus = 360 - eta_plus  % deg, with new diagram
V_plus = sqrt(V_Titan^2 + V_infTitan^2 - ...
              2*V_Titan*V_infTitan*cosd(eta_plus))  % km/s

% Part b

fprintf('\nPart b:\n\n');

xi_plus = asind(V_infTitan/V_plus * sind(eta_plus))  % deg
gamma_plus = xi_plus + gamma_Titan  % deg

h_plus = r*V_plus*cosd(gamma_plus)  % km^2/s
p_plus = h_plus^2 / mu_Saturn  % km
curlE_plus = V_plus^2 / 2 - mu_Saturn / r  % km^2/s^2
a_plus = -mu_Saturn / (2 * curlE_plus)  % km
e_plus = sqrt(1 - p_plus/a_plus)
thetaStar_plus = acosd(1/e_plus * (p_plus/r - 1))  % deg
r_p_plus = a_plus * (1 - e_plus)  % km
r_a_plus = a_plus * (1 + e_plus)  % km
IP_plus = 2*pi*sqrt(a_plus^3 / mu_Saturn)  % s
IP_plus_day = IP_plus/60/60/24  % days
delta_aop = thetaStar_minus - thetaStar_plus  % deg

% Part c

ttl = 'Jordan Mayer - HW 09 Problem 2c';
figure();
plotOrbit2D(e_minus, p_minus, [0 360], 0, '');
hold on;
plotOrbit2D(e_plus, p_plus, [0 360], delta_aop, ttl);
legend('Before Titan encounter', 'After Titan encounter');
axis square;
axis([-3.5 1.5 -3.5 1.5] * 10^6);

% Part d

r_Iapetus = 3561000;  % km
r_Hyperion = 1500000;  % km
ttl = 'Jordan Mayer - HW 09 Problem 2d';
plotOrbit2D(0, r_Iapetus, [0 360], 0, '');
plotOrbit2D(0, r_Hyperion, [0 360], 0, ttl);
legend('s/c before Titan encounter', 's/c after Titan encounter', ...
       'Iapetus', 'Hyperion');
axis([-3.5 3 -3.5 3] * 10^6);

%%% Problem 3 %%

fprintf('\nProblem 3:\n');

% Part a

fprintf('\nPat a:\n\n');

t_Earth = 2453221;  % JD
t_Venus = 2454033;  % JD
t_Mercury = 2454480;  % JD
t_Mercury_i = 2455639;  % JD

TOF_1 = t_Venus - t_Earth  % days
TOF_1_yr = TOF_1/365  % years
TOF_2 = t_Mercury - t_Earth  % days
TOF_2_yr = TOF_2/365  % years
TOF_3 = t_Mercury_i - t_Earth  % days
TOF_3_yr = TOF_3/365  % years

% Part b

fprintf('\nPart b:\n\n');

r_p_minus = r_Mercury  % km
r_a_minus = r_Earth  % km
r = r_Venus  % km
a_minus = 1/2 * (r_p_minus + r_a_minus)  % km
e_minus = 1 - r_p_minus/a_minus
p_minus = a_minus * (1 - e_minus^2)  % km
thetaStar_minus = -acosd(1/e_minus * (p_minus/r - 1))  % deg
E_minus = 2 * atand(((1 - e_minus)/(1 + e_minus))^(1/2) * ...
                    tand(thetaStar_minus/2))  % deg
n_minus = sqrt(mu_Sun/a_minus^3)  % rad/s
M_minus = deg2rad(E_minus) - e_minus * sind(E_minus)
IP_minus = 2*pi/n_minus  % s
IP_minus_day = IP_minus/60/60/24  % days
tMinusTp = M_minus/n_minus + IP_minus  % s
tMinusTp_day = tMinusTp/60/60/24  % days
TOF = tMinusTp_day - 1/2 * IP_minus_day  % days
actual_ratio = (R_Venus + 2990)/R_Venus

% Part c

fprintf('\nPart c:\n\n');

V_Venus = sqrt(mu_Sun/r_Venus)  % km/s

fprintf('\n');

curlE_minus = -mu_Sun/(2*a_minus)  % km^2/s^2
V_minus = sqrt(2 * (curlE_minus + mu_Sun/r))  % km/s
h_minus = sqrt(mu_Sun * p_minus)  % km^2/s
abs_gamma_minus = acosd(h_minus/(r*V_minus))  % deg

fprintf('\n');

V_infVenus = sqrt(V_Venus^2 + V_minus^2 - ...
                  2*V_Venus*V_minus*cosd(abs_gamma_minus))  % km/s
eta_minus = asind(V_Venus/V_infVenus * sind(abs_gamma_minus))  % deg
psi_minus = 180 - abs_gamma_minus - eta_minus  % deg

fprintf('\n');

abs_a = mu_Venus/V_infVenus^2  % km
r_p = 4*R_Venus  % km
e = r_p/abs_a + 1
del = 2*asind(1/e)  % deg

fprintf('\n');

zeta = 1/2 * (180 - del)  % deg
deltaV_eq = sind(del)/sind(zeta) * V_infVenus  % km/s

fprintf('\n');

psi_plus = psi_minus - del  % deg
V_plus = sqrt(V_Venus^2 + V_infVenus^2 - ...
              2 * V_Venus * V_infVenus * cosd(psi_plus))  % km/s
abs_gamma_plus = asind(V_infVenus/V_plus * sind(psi_plus))  % deg

fprintf('\n');

gamma_minus = -abs_gamma_minus  % deg
gamma_plus = -abs_gamma_plus  % deg
h_plus = r*V_plus*cosd(gamma_plus)  % km^2/s
p_plus = h_plus^2/mu_Sun  % km
curlE_plus = V_plus^2 / 2 - mu_Sun/r  % km^2/s^2
a_plus = -mu_Sun/(2*curlE_plus)  % km
e_plus = sqrt(1 - p_plus/a_plus)
thetaStar_plus = -acosd(1/e_plus * (p_plus/r - 1))  % deg
r_p_plus = a_plus * (1 - e_plus)  % km
r_a_plus = a_plus * (1 + e_plus)  % km
IP_plus = 2*pi*sqrt(a_plus^3 / mu_Sun)  % s
IP_plus_day = IP_plus/60/60/24  % days
delta_aop = thetaStar_minus - thetaStar_plus  % deg

% Part d

fprintf('\nPart d:\n\n');

alpha = eta_minus - zeta  % deg

% Part e

fprintf('\nPart e:\n\n');

ttl = 'Jordan Mayer - HW 09, Problem 3d';
figure();
plotOrbit2D(e_minus, p_minus, [0 360], 0, '');
hold on;
plotOrbit2D(e_plus, p_plus, [0 360], delta_aop, '');
plotOrbit2D(0, r_Mercury, [0 360], 0, ttl);
axis square;
axis([-1.5 1.5 -1.5 1.5] * 10^8);
legend('s/c before Venus encounter', 's/c after Venus encounter', ...
       'Mercury');
   
E_plus = 2*atand(((1 - e_plus)/(1 + e_plus))^(1/2) * ...
                 tand(thetaStar_plus/2))  % deg
n_plus = sqrt(mu_Sun/a_plus^3)  % rad/s
M_plus = deg2rad(E_plus) - e_plus * sind(E_plus)
tMinusTp_plus = M_plus/n_plus  % s

fprintf('\n');

thetaStar_1 = -acosd(1/e_plus * (p_plus/r_Mercury - 1))  % deg
E_1 = 2*atand(((1 - e_plus)/(1 + e_plus))^(1/2) * ...
              tand(thetaStar_1/2))  % deg
M_1 = deg2rad(E_1) - e_plus * sind(E_1)
tMinusTp_1 = M_1/n_plus  % s

fprintf('\n');

TOF = tMinusTp_1 - tMinusTp_plus  % s
TOF_day = TOF/60/60/24  % days

% Part f

fprintf('\nPart f:\n\n');

n_Venus = sqrt(mu_Sun/r_Venus^3)  % rad/s
TOF_EarthVenus = 64.2900 * 60 * 60 * 24  % s
TA_EarthVenus = thetaStar_minus + 180  % deg
TA_EarthVenus_rad = deg2rad(TA_EarthVenus)
phi_Venus = TA_EarthVenus_rad - n_Venus * TOF_EarthVenus  % rad
phi_Venus_deg = rad2deg(phi_Venus)  % deg

fprintf('\n');

n_Mercury = sqrt(mu_Sun/r_Mercury^3)  % rad/s
TOF_EarthMercury = TOF_EarthVenus + TOF  % s
TA_EarthMercury = thetaStar_1 + 180  % deg
TA_EarthMercury_rad = deg2rad(TA_EarthMercury)  % rad
phi_Mercury = TA_EarthMercury_rad - n_Mercury * TOF_EarthMercury  % rad
phi_Mercury_deg = rad2deg(phi_Mercury)  % deg