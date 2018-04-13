%%%%%
% AAE 340
% HW 10 - Problem 2
% Jordan Mayer
% 
% Numerically integrate rotational parameters using Euler's EOM in
% body-fixed coordinates for principal axes.
%
% Plot:
%   a) w_x, w_y, w_z (rad/s) vs t
%   b) w_y vs w_x (rad/s)
%   c) H_x, H_y, H_z (kg-m^2/s) vs t
%   d) H_y vs H_z (kg-m^2/s)
%%%%%

close all;   % close any open figures

% Initialize known values
It = 2887;      % principal moment of inertia about transpose axes
                % (x-axis, y-axis), kg-m^2
Ia = 5106;      % principal moment of inertia about axial axis
                % (z-axis), kg-m^2
w_x_0 = 0.6;    % initial angular velocity about x-axis, rad/s
w_y_0 = 0;      % initial angular velocity about y-axis, rad/s
w_z_0 = 1.1;    % initial angular velocity about z-axis, rad/s
lambda = (Ia-It)/It * w_z_0;
wdot_x0 = lambda*w_y_0;    
    % initial change in angular velocity about x-axis, rad/s^2
wdot_y0 = lambda*w_x_0;    
    % initial change in angular velocity about y-axis, rad/s^2
wdot_z0 = 0;    
    % initial change in angular velocity about z-axis, rad/s^2
t0 = 0;         % initial time, s
tf = 10;        % final time, s

% Initialize integration properties
dt = 0.01;      % time step, s
tnext = t0+dt;  % next time, s
options = odeset('RelTol', 1D-12, 'AbsTol', 1E-12, 'InitialStep', dt,...
    'MaxStep', dt);
c = 1;          % counter

% Initialize vectors for state
limit = tf/dt;
state = zeros(limit, 6);    % matrix to hold state variables
y0 = [w_x_0 w_y_0 w_z_0 wdot_x0 wdot_y0 wdot_z0].';     % initial state vector
state(1,:) = y0;

% Initialize vectors for plots
w_x = w_x_0;    % angular velocity about x-axis, rad/s
w_y = w_y_0;    % angular velocity about y-axis, rad/s
w_z = w_z_0;    % angular velocity about z-axis, rad/s
t = t0:dt:tf;   % time, s

% Numerically integrate for system
tic;    % start timer to check for excessively long loop
while tnext <= tf
    c = c + 1;    % increment counter
    [tn, yn] = ode45('Euler_PMOI', [t0 tnext], y0, options);    % integrate
    % set up next
    t0 = tn(end);
    tnext = tnext + dt;
    y0 = yn(end, 1:6);
    state(c,:) = y0;    % save state at this time
    w_x(c) = y0(1);
    w_y(c) = y0(2);
    w_z(c) = y0(3);
    if (toc > 60*5)   % throw error if in loop longer than 5 minutes
        error('We''ve been looping for 5 minutes. Speed it up!');
    end
end

% Calculate angular momentum
Hx = It * w_x;  % angular momentum about x-axis, kg-m^2/s
Hy = It * w_y;  % angular momentum about y-axis, kg-m^2/s
Hz = Ia * w_z;  % angular momentum about z-axis, kg-m^2/s

% Create and format plots
% plot a
figure('Units', 'inches', 'Position', [1 4 5 4]); hold on;
plot(t, w_x, '-r'); plot(t, w_y, '-g'); plot(t, w_z, '-b');
title('Problem 2a'); xlabel('time, s'); ylabel('angular velocity, rad/s');
legend('w_x', 'w_y', 'w_z'); grid on; hold off;
% plot b
figure('Units', 'inches', 'Position', [1 4 5 4]); plot(w_x, w_y, '-k');
title('Problem 2b: angular velocity');
xlabel('w_x, rad/s'); ylabel('w_y, rad/s'); grid on;
% plot c
figure('Units', 'inches', 'Position', [1 4 5 4]); hold on;
plot(t, Hx, '-r'); plot(t, Hy, '-g'); plot(t, Hz, '-b');
title('Problem 2c'); xlabel('time, s'); 
ylabel('angular momentum, kg-m^2/s');
legend('H_x', 'H_y', 'H_z'); grid on; hold off;
% plot d
figure('Units', 'inches', 'Position', [1 4 5 4]); plot(Hx, Hy, '-k');
title('Problem 2d: angular momentum');
xlabel('H_x, rad/s'); ylabel('H_y, rad/s'); grid on;