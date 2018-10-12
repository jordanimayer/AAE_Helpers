%%%%%
% Jordan Mayer
% AAE 532
% HW 05
%%%%%

%%% Problem 2 %%%

fprintf('\nProblem 2:\n\n');

% Part a
fprintf('Part a:\n\n');
% Find E from M
M_0 = deg2rad(360-120);  % initial mean anomaly
E_0 = 123;  % indicates E_0 was not found
e = 0.632653;
for E = linspace(0,2*pi,10^6)  % iterate through possible eccentric
                               % anomaly values
    M = E - e*sin(E);
    if abs(M - M_0) < 10^-5  % stop when value of E results in (almost)
                             % correct M
        E_0 = E;
        break;
    end
end
fprintf('E_0 = %.6f\n\n', E_0);

% Part c
fprintf('Part b:\n\n');
% Find E from M
M_f = 2.61798;  % final mean anomaly
E_f = 123;  % indicates E_f was not found
for E = linspace(0,2*pi,10^6)  % iterate through possible eccentric
                               % anomaly values
    M = E - e*sin(E);
    if abs(M - M_f) < 10^-5  % stop when value of E results in (almost)
                             % correct M
        E_f = E;
        break;
    end
end
fprintf('E_f = %.6f\n\n', E_f);

% Part d
p = 9980.64;  % semi-latus rectum, km
thetaStar = linspace(rad2deg(-2.81970),rad2deg(2.98740),1000);  % range of true anomalies, deg
r = p ./ (1 + e*cosd(thetaStar));  % corresponding radii, km
r_1 = r .* cosd(thetaStar);  % radii components in e_hat dir
r_2 = r .* sind(thetaStar);  % radii components in p_hat dir

close all;  % close any open figures
plot(r_1, r_2); grid on;
xlabel('r_1 (km)'); ylabel('r_2 (km)');
title('Jordan Mayer - HW 05 Problem 2 Part d');
axis([-30000 10000 -20000 20000]);  % scale axes for more intuitive plot