%%%%%
% Jordan Mayer
% AAE 532
% HW 10
%
% Perform arithmetic to assist in orbital analysis
%%%%%

%%% Preliminary setup %%%

clear all; close all; format compact; format long; rehash;

global R_Earth mu_Earth r_Earth ...
       R_Mars mu_Mars r_Mars ...
       mu_Sun;
  % descriptions and units in load_constants.m

load_constants;

%%% Problem 1 %%%

fprintf('\nProblem 1:\n');

% Part a

fprintf('\nPart a:\n\n');

r_1 = 2.5*R_Earth  % km
V_1_minus = sqrt(mu_Earth/r_1)  % km/s
r_2 = 6*R_Earth  % km
a_T = 1/2 * (r_1 + r_2)  % km
curlE_T = -mu_Earth/(2*a_T)  % km^2/s^2
V_1_plus = sqrt(2*(curlE_T + mu_Earth/r_1))  % km/s
deltaV_1 = V_1_plus - V_1_minus  % km/s

fprintf('\n');

V_2_minus = sqrt(2*(curlE_T + mu_Earth/r_2))  % km/s
V_2_plus = sqrt(mu_Earth/r_2)  % km/s
deltaV_2 = V_2_plus - V_2_minus  % km/s

fprintf('\n');

deltaV_Tot = deltaV_2 + deltaV_2  % km/s

fprintf('\n');

TOF = pi * sqrt(a_T^3 / mu_Earth)  % s
TOF_hr = TOF/60/60  % hours

fprintf('\n');

n_2 = sqrt(mu_Earth / r_2^3)  % rad/s
phi = pi - n_2 * TOF  % rad
phi_deg = rad2deg(phi)

fprintf('\n');

n_1 = sqrt(mu_Earth / r_1^3)  % rad/s
t_s = 2*pi/(n_1 - n_2)  % s
t_s_hr = t_s/60/60  % hours

% Part b

fprintf('\nPart b:\n\n');

c = sqrt(r_1^2 + r_2^2)  % km
s = (r_1 + r_2 + c)/2  % km
a = s/2  % km
alpha = 2*asind(sqrt(s/(2*a)))  % deg
beta = 2*asind(sqrt((s-c)/(2*a)))  % deg
p = 4*a*(s-r_1)*(s-r_2)/c^2 * (sind((alpha + beta)/2))^2  % km
p = 4*a*(s-r_1)*(s-r_2)/c^2 * (sind((alpha - beta)/2))^2 % km
  % should be the same (it is!)
e = sqrt(1 - p/a)
r_p = a*(1-e)  % km

fprintf('\n');

TOF = sqrt(a^3/mu_Earth) * ((deg2rad(alpha) - deg2rad(beta)) - ...
      (sind(alpha) - sind(beta)))  % s
TOF_hr = TOF/60/60  % hours
r_A = r_2
r_D = r_1
curlE = -mu_Earth/(2*a)  % km^2/s^2
V_D = sqrt(2*(curlE + mu_Earth/r_D))  % km/s
V_A = sqrt(2*(curlE + mu_Earth/r_A))  % km/s
thetaStar_D = acosd(1/e * (p/r_D - 1))  % deg
thetaStar_A = acosd(1/e * (p/r_A - 1))  % deg
h = sqrt(mu_Earth * p)  % km^2/s
gamma_D = -acosd(h/(r_D * V_D))  % deg
gamma_A = acosd(h/(r_A * V_A))  % deg

fprintf('\n');

TA = deg2rad(270)  % rad
phi = TA - n_2 * TOF  % rad
phi_deg = rad2deg(phi)  % deg

% Part c

fprintf('\nPart c:\n\n');

V_1 = V_1_minus  %  km, from Part a
deltaV_D = sqrt(V_1^2 + V_D^2 - 2*V_1*V_D*cosd(gamma_D))  % km/s
beta_1 = asind(V_D/deltaV_D * sind(abs(gamma_D)))  % deg
beta_2 = 180 - beta_1  % deg
eta_1 = asind(V_1/deltaV_D * sind(abs(gamma_D)))  % deg
eta_2 = 180 - eta_1  % deg
beta = beta_1  % deg
eta = eta_1  % deg
abs_alpha_D = 180 - beta  % deg
alpha_D = -abs_alpha_D

fprintf('\n');

deltaV_D_V = deltaV_D*cosd(alpha_D)  % km/s
deltaV_D_C = deltaV_D*sind(alpha_D)  % km/s

fprintf('\n');

V_2 = V_2_plus  % km, from Part a
deltaV_A = sqrt(V_2^2 + V_A^2 - 2*V_2*V_A*cosd(gamma_A))  % km/s
beta_1 = asind(V_2/deltaV_A * sind(gamma_A))  % deg
beta_2 = 180 - beta_1  % deg
eta_1 = asind(V_A/deltaV_A * sind(gamma_A))  % deg
eta_2 = 180 - eta_1  % deg
beta = beta_2  % deg
eta = eta_2  % deg
abs_alpha_A = 180 - beta  % deg
alpha_A = -abs_alpha_A  % deg

fprintf('\n');

deltaV_A_V = deltaV_A*cosd(alpha_A)  % km/s
deltaV_A_C = deltaV_A*sind(alpha_A)  % km/s

fprintf('\n');

deltaV_Tot = deltaV_D + deltaV_A  % km/s, for part f

% Part d

fprintf('\nPart d:\n\n');

IP_1 = 2*pi*sqrt(r_1^3/mu_Earth)  % s
IP_2 = 2*pi*sqrt(r_2^3/mu_Earth)  % s
t_f = 2 * (IP_1 + IP_2) + TOF

% Part e

fprintf('\nPart e:\n');

% Part e-b

fprintf('\nPart e-b:\n\n');

TA = 150;  % deg
c = sqrt(r_1^2 + r_2^2 - 2*r_1*r_2*cosd(TA))  % km
s = (r_1 + r_2 + c)/2  % km
a = s/2  % km
alpha = 2*asind(sqrt(s/(2*a)))  % deg
beta = 2*asind(sqrt((s-c)/(2*a)))  % deg
p = 4*a*(s-r_1)*(s-r_2)/c^2 * (sind((alpha + beta)/2))^2  % km
e = sqrt(1 - p/a)
r_p = a*(1-e)  % km

fprintf('\n');

TOF = sqrt(a^3/mu_Earth) * ((deg2rad(alpha) - deg2rad(beta)) - ...
      (sind(alpha) - sind(beta)))  % s
TOF_hr = TOF/60/60  % hours
r_D = r_1  % km
r_A = r_2  % km
thetaStar_D = acosd(1/e * (p/r_D - 1))  % deg
thetaStar_A = acosd(1/e * (p/r_A - 1))  % deg
curlE = -mu_Earth/(2*a)  % km^2/s^2
V_D = sqrt(2*(curlE + mu_Earth/r_D))  % km/s
V_A = sqrt(2*(curlE + mu_Earth/r_A))  % km/s
h = sqrt(mu_Earth * p)  % km^2/s
gamma_D = acosd(h/(r_D * V_D))  % deg
gamma_A = -acosd(h/(r_A * V_A))  % deg

fprintf('\n');

TA = deg2rad(150)  % rad
phi = TA - n_2 * TOF  % rad
phi_deg = rad2deg(phi)  % deg

% Part e-c

fprintf('\nPart e-c:\n\n')

deltaV_D = sqrt(V_1^2 + V_D^2 - 2*V_1*V_D*cosd(abs(gamma_D)))  % km/s
beta1 = asind(V_D/deltaV_D * sind(abs(gamma_D)))  % deg
beta2 = 180 - beta1
eta1 = asind(V_1/deltaV_D * sind(abs(gamma_D)))  % deg
eta2 = 180 - eta1
beta = beta2  % deg
eta = eta1  % deg
abs_alpha_D = 180 - beta  % deg
alpha_D = abs_alpha_D  % deg

fprintf('\n')

deltaV_D_V = deltaV_D * cosd(alpha_D)  % km/s
deltaV_D_C = deltaV_D * sind(alpha_D)  % km/s

fprintf('\n')

deltaV_A = sqrt(V_A^2 + V_2^2 - 2*V_A*V_2*cosd(abs(gamma_A)))  % km/s
beta1 = asind(V_2/deltaV_A * sind(abs(gamma_A)))  % deg
beta2 = 180 - beta1  % deg
eta1 = asind(V_A/deltaV_A * sind(abs(gamma_A)))  % deg
eta2 = 180 - eta1  % deg
beta = beta2  % deg
eta = eta1  % deg
abs_alpha_A = 180 - beta  % deg
alpha_A = abs_alpha_A  % deg

fprintf('\n');

deltaV_A_V = deltaV_A*cosd(alpha_A)  % km/s
deltaV_A_C = deltaV_A*sind(alpha_A)  % km/s

fprintf('\n');

deltaV_Tot = deltaV_D + deltaV_A  % km/s, for part f

% Part e-d

fprintf('\nPart e-d:\n\n');

t_f = 2*(IP_1 + IP_2) + TOF  % s

%%% Problem 2 %%%

fprintf('\nProblem 2:\n');

% Part a-i

fprintf('\nPart a-i:\n\n');

a = 1/2 * (r_Earth + r_Mars)  % km
TOF_Hoh = pi * sqrt(a^3/mu_Sun)  % s
TOF_Hoh_day = TOF/60/60/24  % days

fprintf('\n');
n_Earth = sqrt(mu_Sun/r_Earth^3)  % rad/s
n_Mars = sqrt(mu_Sun/r_Mars^3)  % rad/s
t_s = 2*pi/(n_Earth - n_Mars)  % s
t_s_day = t_s/60/60/24  % days

% Part a-ii

fprintf('\nPart a-ii:\n\n');

V_Earth = sqrt(mu_Sun/r_Earth)  % km/s
curlE = -mu_Sun/(2*a)  % km^2/s^2
V_plus = sqrt(2*(curlE + mu_Sun/r_Earth))  % km/s
V_infEarth = V_plus - V_Earth  % km/s

fprintf('\n');

r_i = 250 + R_Earth  % km
V_1_minus = sqrt(mu_Earth/r_i)  % km/s
curlE_1_plus = V_infEarth^2 / 2  % km^2/s^2
V_1_plus = sqrt(2*(curlE_1_plus + mu_Earth/r_i))  % km/s
deltaV_1 = V_1_plus - V_1_minus  % km/s

fprintf('\n');

V_Mars = sqrt(mu_Sun/r_Mars)  % km/s
V_minus = sqrt(2*(curlE + mu_Sun/r_Mars))  % km/s
V_infMars = V_Mars - V_minus  % km/s

fprintf('\n');

r_f = 1000 + R_Mars  % km
V_2_plus = sqrt(mu_Mars/r_f)  % km/s
curlE_2_minus = V_infMars^2 / 2  % km^2/s^2
V_2_minus = sqrt(2*(curlE_2_minus + mu_Mars/r_f))  % km/s
deltaV_2 = V_2_minus - V_2_plus  % km/s

fprintf('\n');

deltaV_Tot = deltaV_1 + deltaV_2  % km/s

% Part b

fprintf('\nPart b:\n\n');

TA = 160;  % deg
r_1 = r_Earth  % km
r_2 = r_Mars  % km
c = sqrt(r_1^2 + r_2^2 - 2*r_1*r_2*cosd(TA))  % km
s = (r_1 + r_2 + c)/2  % km
a = s/2  % km
alpha = 2*asind(sqrt(s/(2*a)))  % deg
beta = 2*asind(sqrt((s-c)/(2*a)))  % deg
beta_alt = 2*(180 - asind(sqrt((s-c)/(2*a))))  % deg
p = 4*a*(s-r_1)*(s-r_2)/c^2 * (sind((alpha + beta)/2))^2  % km
p = 4*a*(s-r_1)*(s-r_2)/c^2 * (sind((alpha - beta)/2))^2  % km
p = 4*a*(s-r_1)*(s-r_2)/c^2 * (sind((alpha + beta_alt)/2))^2  % km
p = 4*a*(s-r_1)*(s-r_2)/c^2 * (sind((alpha - beta_alt)/2))^2  % km
  % should all be the same (they are!)
e = sqrt(1 - p/a)
r_p = a*(1-e)  % km
r_a = a*(1+e)  % km
curlE = -mu_Sun/(2*a)  % km^2/s^2

fprintf('\n');

thetaStar_D = acosd(1/e * (p/r_1 - 1))  % deg
thetaStar_A = -acosd(1/e * (p/r_2 - 1))  % deg
V_D = sqrt(2*(curlE + mu_Sun/r_1))  % km/s
V_A = sqrt(2*(curlE + mu_Sun/r_2))  % km/s
h = sqrt(mu_Sun * p)  % km^2/s
gamma_D = acosd(h/(r_1*V_D))  % deg
gamma_A = -acosd(h/(r_2*V_A))  % deg
TOF_min = sqrt(a^3 / mu_Sun) * ((deg2rad(alpha) - sind(alpha)) - ...
                                (deg2rad(beta) - sind(beta)))  % s
TOF_min_day = TOF_min/60/60/24  % days

% plot initial orbit
plotOrbit2D(0, r_1, [0 360], 0, '');
hold on;
% plot transfer orbit
plotOrbit2D(e, p, [thetaStar_D, 360 + thetaStar_A], 0, '');
% plot final orbit
ttl = 'Jordan Mayer - HW 10, Problem 2b';
plotOrbit2D(0, r_2, [0 360], 0, ttl);
legend('Earth Orbit', 'Transfer', 'Mars Orbit');
axis square; axis([-2.5 2.5 -2.5 2.5] * 10^8);

% Part c

fprintf('\nPart c:\n\n');

phi = deg2rad(TA) - n_Mars * TOF_min  % rad
phi_deg = rad2deg(phi)  % deg

% Part d

fprintf('\nPart d:\n\n');

V_infEarth = sqrt(V_Earth^2 + V_D^2 - ...
                  2*V_Earth*V_D*cos(abs(gamma_D)))  % km/s
V_1_minus = V_1_minus  % km/s
curlE_1_plus = V_infEarth^2 / 2  % km^2/s^2
V_1_plus = sqrt(2*(curlE_1_plus + mu_Earth/r_i))  % km/s
deltaV_D = V_1_plus - V_1_minus  % km/s

fprintf('\n');

V_infMars = sqrt(V_Mars^2 + V_A^2 - ...
                 2*V_Mars*V_A*cosd(abs(gamma_A)))  % km/s

fprintf('\n');

V_2_plus = V_2_plus  % km/s
curlE_2_minus = V_infMars^2 / 2  % km^2/s^2
V_2_minus = sqrt(2*(curlE_2_minus + mu_Mars/r_f))  % km/s
deltaV_A = V_2_minus - V_2_plus  % km/s

% Part e

fprintf('\nPart e:\n\n');

deltaV_Tot = deltaV_D + deltaV_A  % km/s