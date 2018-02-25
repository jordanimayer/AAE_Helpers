% AAE 333: HW 3 - Problem 2c
%
% This script lots density and pressure as functions of altitude from 0 to
% 10 km
%
% @author Jordan Mayer, mayer15
% @version 2.1.2017

h = 0:0.05:10;    %altitude, km
T = 15 - 65*h/ 10;       %temperature, degress Celsius
p0 = 1.01 * 10^5;       %pressure at sea-level, Pa
R = 287.1;              %ideal gas constant for air, J/(kg*K)
g = 9.81;               %gravitational constant, m/s^2

rho = p0 ./ (R * (T + 273) + g * (h * 1000));  %density, kg/m^3
p = rho * R .* (T + 273);       %pressure, Pa

%plot density vs altitude
figure(1);
plot(rho, h);
xlabel('Density, kg/m^3');
ylabel('Altitude, km');
title('Standard Atmosphere Approximation');
grid on;

%plot pressure vs altitude
figure(2);
plot(p, h);
xlabel('Pressure, Pa');
ylabel('Altitude, km');
title('Standard Atmosphere Approximation');
grid on;

stand_atmo = [h(1) T(1) p(1)*10^(-4) rho(1); h(21) T(21) p(21)*10^(-4) rho(21); h(41) T(41) p(41)*10^(-4) rho(41); h(61) T(61) p(61)*10^(-4) rho(61); h(81) T(81) p(81)*10^(-4) rho(81); h(101) T(101) p(101)*10^(-4) rho(101); h(121) T(121) p(121)*10^(-4) rho(121);h(141) T(141) p(141)*10^(-4) rho(141); h(161) T(161) p(161)*10^(-4) rho(161); h(181) T(181) p(181)*10^(-4) rho(181);h(201) T(201) p(201)*10^(-4) rho(201)]
