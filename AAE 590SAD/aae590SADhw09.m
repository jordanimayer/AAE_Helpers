%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 09
%
% Perform numerical simulations and plot results.
%%%%%

%% Preliminary setup

close all; clear all; format compact; format short; rehash toolbox;
addpath('../helpers');

% Problem 2

% numerically simulate dynamic and kinematic differential equations
% simultaneously

%% Part ci
close all; clear all;

% set up constants
sigma_0 = 6;  % deg
v_f = 4;  % final number of revs
ttl_start = 'Jordan Mayer, AAE 590/440, HW 09, ';

% plot particular solution for nutation angle (all zero)
sigma_part = [0, 0];
figure(1);
plot([0, v_f], sigma_part); hold on;

for r = [1.8, 22]
    hw09_sim_and_plot(v_f, sigma_0, r, 1, 1);
end

% format plots
figure(1);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem ci: Numerical Results'];
title(ttl); grid on;
legend('Particular Solution', 'r = 1.8', 'r = 22');
figure(2);
xlabel('Number of Revolutions'); ylabel('Error');
title(ttl); grid on;
legend('r = 1.8', 'r = 22');

%% Part cii
close all; clear all;

% set up constants
v_f = 4;  % final number of revs
options = odeset('RelTol', 1.0e-7, 'AbsTol', 1.0e-7);
  % set error tolerance to 10^(-7)
ttl_start = 'Jordan Mayer, AAE 590/440, HW 09, ';

% plot particular solution for nutation angle (all zero)
sigma_part = [0, 0];
figure(3);
plot([0, v_f], sigma_part); hold on;
figure(5);
plot([0, v_f], sigma_part); hold on;
figure(7);
plot([0, v_f], sigma_part); hold on;

for sigma_0 = [3, 6]  % deg
    for r = [1.5, 8, 22]
        if r == 1.5
            fig = 3;
        elseif r == 8
            fig = 5;
        else
            fig = 7;
        end
        
        hw09_sim_and_plot(v_f, sigma_0, r, fig, 1);
    end
end

% format plots
% r = 1.5
figure(3);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem cii: Numerical Results, r = 1.5'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');
figure(4);
xlabel('Number of Revolutions'); ylabel('Error');
title(ttl); grid on;
legend('sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = 8
figure(5);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem cii: Numerical Results, r = 8'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');
figure(6);
xlabel('Number of Revolutions'); ylabel('Error');
title(ttl); grid on;
legend('sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = 22
figure(7);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem cii: Numerical Results, r = 22'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');
figure(8);
xlabel('Number of Revolutions'); ylabel('Error');
title(ttl); grid on;
legend('sigma_0 = 3 deg', 'sigma_0 = 6 deg');

%% Part ciii

% set up constants
v_0 = 0;  % initial number of revs
options = odeset('RelTol', 1.0e-7, 'AbsTol', 1.0e-7);
  % set error tolerance to 10^(-7)
ttl_start = 'Jordan Mayer, AAE 590/440, HW 09, ';

% copy cii plots
figure(9);
copyobj(allchild(3), 9);
figure(10);
copyobj(allchild(5), 10);
figure(11);
copyobj(allchild(7), 11);

sigma_0 = 12;  % deg

for r = [1.5, 8, 22]
    if r == 1.5
        fig = 9;
    elseif r == 8
        fig = 10;
    else
        fig = 11;
    end
    
    hw09_sim_and_plot(v_f, sigma_0, r, fig, 0);
end

% format plots
figure(9);
title([ttl_start, 'Problem ciii: Numerical Results, r = 1.5']);
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg', ...
       'sigma_0 = 12 deg');

% r = 8
figure(10);
title([ttl_start, 'Problem ciii: Numerical Results, r = 8']);
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg', ...
       'sigma_0 = 12 deg');

% r = 22
figure(11);
title([ttl_start, 'Problem ciii: Numerical Results, r = 22']);
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg', ...
       'sigma_0 = 12 deg');
   
%% Part d

%% try more simulation time (20 revs instead of 4)
close all; clear all;

% set up constants
v_f = 20;  % final number of revs
options = odeset('RelTol', 1.0e-7, 'AbsTol', 1.0e-7);
  % set error tolerance to 10^(-7)
ttl_start = 'Jordan Mayer, AAE 590/440, HW 09, ';

% plot particular solution for nutation angle (all zero)
sigma_part = [0, 0];
figure(1);
plot([0, v_f], sigma_part); hold on;
figure(2);
plot([0, v_f], sigma_part); hold on;
figure(3);
plot([0, v_f], sigma_part); hold on;

for sigma_0 = [3, 6]  % deg
    for r = [1.5, 8, 22]
        if r == 1.5
            fig = 1;
        elseif r == 8
            fig = 2;
        else
            fig = 3;
        end
        
        hw09_sim_and_plot(v_f, sigma_0, r, fig, 0);
    end
end

% format plots
% r = 1.5
figure(1);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = 1.5'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = 8
figure(2);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = 8'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = 22
figure(3);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = 22'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

%% try new r values
close all; clear all;

% set up constants
v_f = 20;  % final number of revs
options = odeset('RelTol', 1.0e-7, 'AbsTol', 1.0e-7);
  % set error tolerance to 10^(-7)
ttl_start = 'Jordan Mayer, AAE 590/440, HW 09, ';

% plot particular slution for nutation angle (all zero)
sigma_part = [0, 0];
figure(1);
plot([0, v_f], sigma_part); hold on;
figure(2);
plot([0, v_f], sigma_part); hold on;

for sigma_0 = [3, 6]  % deg
    for r = [-1.8, 36]
        if r == -1.8
            fig = 1;
        else
            fig = 2;
        end
        
        hw09_sim_and_plot(v_f, sigma_0, r, fig, 0);
    end
end

% format plots
% r = -1.8
figure(1);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = -1.8'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = 36
figure(2);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = 36'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

%% check symmetry
close all; clear all;

% set up constants
v_f = 20;  % final number of revs
options = odeset('RelTol', 1.0e-7, 'AbsTol', 1.0e-7);
  % set error tolerance to 10^(-7)
ttl_start = 'Jordan Mayer, AAE 590/440, HW 09, ';

% plot particular slution for nutation angle (all zero)
sigma_part = [0, 0];
figure(1);
plot([0, v_f], sigma_part); hold on;
figure(2);
plot([0, v_f], sigma_part); hold on;
figure(3);
plot([0, v_f], sigma_part); hold on;
figure(4);
plot([0, v_f], sigma_part); hold on;

for sigma_0 = [3, 6]  % deg
    for r = [2, -2, 20, -20]
        if r == 2
            fig = 1;
        elseif r == -2
            fig = 2;
        elseif r == 20
            fig = 3;
        else
            fig = 4;
        end
        
        hw09_sim_and_plot(v_f, sigma_0, r, fig, 0);
    end
end

% format plots
% r = +-2
figure(1);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = 2'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = -2
figure(2);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = -2'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = 20
figure(3);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = 20'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');

% r = -20
figure(4);
xlabel('Number of Revolutions'); ylabel('Nutation Angle, deg');
ttl = [ttl_start, 'Problem d: Numerical Results, r = -20'];
title(ttl); grid on;
legend('Particular Solution', 'sigma_0 = 3 deg', 'sigma_0 = 6 deg');