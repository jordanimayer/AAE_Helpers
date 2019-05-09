%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 10
%
% Perform numerical simulations and plot results.
% (Also do some basic arithmetic).
%%%%%

%% Preliminary setup

close all; clear all; format compact; format short; rehash toolbox;
addpath('../helpers');

%% Problem 1i
fprintf('\nProblem 1i\n\n');

mu_Earth = 3.986E5;  % Earth gravitational parameter, km^3/s^2
R_Earth = 6378;  % Earth radius, km
r = 100;

for h = [200, 800, 2000, 25000]  % orbit altitude, km
    fprintf('h =\n');
    disp(h);
    omega = sqrt(mu_Earth/(R_Earth + h)^3)  % rad/s
    angVel_30 = r*omega  % rad/s
    fprintf('\n');
end

%% Problem 1ii
close all; clear all;

% set up constants
sigma_0 = 6;  % deg
v_span = [0, 2];  % span of revolutions to simulate for
r_list = [1.5, 8, 22, 120, 200];
ttl = 'Jordan Mayer, AAE 590/440, HW10, Problem 1, Numerical Results';
lgd = {'particular solution', 'r = 1.5', 'r = 8', 'r = 22', 'r = 120', ...
       'r = 200'};
lgd_lambda = {'particular solution', 'r = 1.5', 'r = 1.5 (slope)', ...
              'r = 8', 'r = 8 (slope)', 'r = 22', 'r = 22 (slope)', ...
              'r = 120', 'r = 120 (slope)'};
xlab = 'Number of Revolutions';
fpos = [200 300 650 500];

% plot particular solutions (sigma = eta = 0, lambda = 2pi*v rads)
angle_part = [0, 0];
f1 = figure(1); plot(v_span, angle_part); hold on;
  % for sigma results
f2 = figure(2); plot(v_span, angle_part); hold on;
  % for eta results
f3 = figure(3); plot(v_span, rad2deg(v_span*2*pi)); hold on;
  % for lambda results
f4 = figure(4); plot(v_span, rad2deg(v_span*2*pi)); hold on;
  % for delta_lambda
set(f1, 'Position', fpos);
set(f2, 'Position', fpos);
set(f3, 'Position', fpos);
% simulate and plot angles
for r = r_list
    hw10_sim_and_plot(v_span, sigma_0, r, 1);
end

% format plots
figure(1);
xlabel(xlab); ylabel('Nutation Angle (sigma), deg');
title(ttl); grid on; legend(lgd);
figure(2);
xlabel(xlab); ylabel('Orbital Precession Angle (eta), deg');
title(ttl); grid on; legend(lgd);
figure(3);
xlabel(xlab); ylabel('Inertial Precession Angle (lambda), deg');
title(ttl); grid on; legend(lgd_lambda);

%% Problem 1vi
close all; clear all;

% set up constants
sigma_0 = 6;  % deg
v_span = [0, 2];  % span of revolutions to simulate for
r_list = -1*[1.5, 8, 22, 120];
ttl = 'Jordan Mayer, AAE 590/440, HW10, Problem 1, Numerical Results';
lgd = {'particular solution', 'r = -1.5', 'r = -8', 'r = -22', ...
       'r = -120'};
xlab = 'Number of Revolutions';
fpos = [200 300 650 500];

% plot particular solutions (sigma = eta = 0, lambda = 2pi*v rads)
angle_part = [0, 0];
f1 = figure(1); plot(v_span, angle_part); hold on;
  % for sigma results
f2 = figure(2); plot(v_span, angle_part); hold on;
  % for eta results
f3 = figure(3); plot(v_span, rad2deg(v_span*2*pi)); hold on;
  % for lambda results
set(f1, 'Position', fpos);
set(f2, 'Position', fpos);
set(f3, 'Position', fpos);

% simulate and plot angles
for r = r_list
    hw10_sim_and_plot(v_span, sigma_0, r, 1);
end

% format plots
figure(1);
xlabel(xlab); ylabel('Nutation Angle (sigma), deg');
title(ttl); grid on; legend(lgd);
figure(2);
xlabel(xlab); ylabel('Orbital Precession Angle (eta), deg');
title(ttl); grid on; legend(lgd);
figure(3);
xlabel(xlab); ylabel('Inertial Precession Angle (lambda), deg');
title(ttl); grid on; legend(lgd);

%% Problem 2d
fprintf('\nProblem 2d:\n\n');
close all; clear all;

% find characteristic roots for different r values
I = 500;  % kg*m^2
J = 150;  % kg*m^2
x = (J-I)/I;
sigma_0 = 6;  % deg
r_list = [-20, -1.5, 0, 1.5, 8, 22, 120];

for r = r_list
    y = r - 1;
    Q = y + x*(y+1);
    b = 3*x + Q^2 + 1;
    c = 3*x*Q + Q^2;
      % for characteristic equation: lambda^4 + b*lambda^2 + c = 0
    lambda_nd = zeros(1,4);
      % characteristic roots, nondimensionalized w.r.t. orbital rate
    
    lambda_nd(1) = sqrt((-b + sqrt(b^2 - 4*c))/2);
    lambda_nd(2) = - lambda_nd(1);
    lambda_nd(3) = sqrt((-b - sqrt(b^2 - 4*c))/2);
    lambda_nd(4) = -lambda_nd(3);
    
    fprintf('Characteristic roots for r = %.1f:\n', r);
    disp(lambda_nd(1));
    disp(lambda_nd(2));
    disp(lambda_nd(3));
    disp(lambda_nd(4));
end

%% get numerical results

% set up constants for numerical integration
v_span = [0, 2];  % span of revolutions to simulate for
sigma0_list = [3, 6];  % deg
ttl = 'Jordan Mayer, AAE 590/440, HW10, Problem 2, Numerical Results, ';
lgd = {'particular solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg'};
xlab = 'Number of Revolutions'; ylab = 'Nutation Angle, deg';
fpos = [200 300 650 500];
sigma_part = [0, 0];

fignum = 0;
for r = r_list
    fignum = fignum + 1;
    
    % plot particular solution
    fg1 = figure(fignum); hold off; plot(v_span, sigma_part); hold on;
    set(fg1, 'Position', fpos);
    
    % plot comparison solutions
    for sigma_0 = sigma0_list
        hw10_sim_and_plot(v_span, sigma_0, r, fignum);
    end
    
    % format plot
    figure(fignum);
    xlabel(xlab); ylabel(ylab);
    title([ttl, 'r = ', num2str(r)]); grid on;
    legend(lgd);
end

%% determine boundary r values

r = 0; pos_stbl = 0; neg_stbl = 0;

while pos_stbl == 0
    r = r + 0.001;
    y = r - 1;
    Q = y + x*(y+1);
    b = 3*x + Q^2 + 1;
    c = 3*x*Q + Q^2;
      % for characteristic equation: lambda^4 + b*lambda^2 + c = 0
    lambda_nd = zeros(1,4);
      % characteristic roots, nondimensionalized w.r.t. orbital rate
    
    lambda_nd(1) = sqrt((-b + sqrt(b^2 - 4*c))/2);
    lambda_nd(2) = - lambda_nd(1);
    lambda_nd(3) = sqrt((-b - sqrt(b^2 - 4*c))/2);
    lambda_nd(4) = -lambda_nd(3);
    
    if size(find(real(lambda_nd) > 0), 2) == 0
        % no lambda_nd has positive real part --> not unstable!
        pos_stbl = 1;
        r_crit_pos = r
    end
    
    if r > 1000
        error('r is greater than 1000 and still unstable!');
    end
end

r = 0;

while neg_stbl == 0
    r = r - 0.001;
    y = r - 1;
    Q = y + x*(y+1);
    b = 3*x + Q^2 + 1;
    c = 3*x*Q + Q^2;
      % for characteristic equation: lambda^4 + b*lambda^2 + c = 0
    lambda_nd = zeros(1,4);
      % characteristic roots, nondimensionalized w.r.t. orbital rate
    
    lambda_nd(1) = sqrt((-b + sqrt(b^2 - 4*c))/2);
    lambda_nd(2) = - lambda_nd(1);
    lambda_nd(3) = sqrt((-b - sqrt(b^2 - 4*c))/2);
    lambda_nd(4) = -lambda_nd(3);
    
    if size(find(real(lambda_nd) > 0), 2) == 0
        % no lambda_nd has positive real part --> not unstable!
        neg_stbl = 1;
        r_crit_neg = r
    end
    
    if r > 1000
        error('r is greater than 1000 and still unstable!');
    end
end

%% Problem 2e
close all; clear all;

% simulate for new r values
I = 500;  % kg*m^2
J = 150;  % kg*m^2
x = (J-I)/I;
sigma_0 = 6;  % deg
r_list = [5, 12];

% set up constants for numerical integration
v_span = [0, 10];  % span of revolutions to simulate for
sigma0_list = [3, 6];  % deg
ttl = 'Jordan Mayer, AAE 590/440, HW10, Problem 2, Numerical Results, ';
lgd = {'particular solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg'};
xlab = 'Number of Revolutions';
fpos = [200 300 650 500];
sigma_part = [0, 0];
lambda_part = 2*pi*v_span;

fignum = -1;
for r = r_list
    fignum = fignum + 2;
    
    % plot particular solution
    fg1 = figure(fignum); hold off; plot(v_span, sigma_part); hold on;
    fg2 = figure(fignum+1); hold off; plot(v_span, lambda_part); hold on;
    set(fg1, 'Position', fpos); set(fg2, 'Position', fpos);
    
    % plot comparison solutions
    for sigma_0 = sigma0_list
        hw10_sim_and_plot(v_span, sigma_0, r, fignum);
    end
    
    % format plots
    figure(fignum);
    xlabel(xlab); ylabel('Nutation Angle, deg');
    title([ttl, 'r = ', num2str(r)]); grid on;
    legend(lgd); ylim([0 11]);
    figure(fignum + 1);
    xlabel(xlab); ylabel('Inertial Precession Angle, deg');
    title([ttl, 'r = ', num2str(r)]); grid on;
end