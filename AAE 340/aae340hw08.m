 %%%%%
% AAE 340
% HW 08
% Jordan Mayer
%%%%%

%%% Problem 1 %%%
% Numerically integrate given system for 3 cycles of motion.
% 
% Part a
% Find time for 3 cycles.
%
% Part b
% Generate the following plots:
%   i. theta (deg) as a function of time
%   ii. x as a function of time
%   iii. theta (deg) as a function of x
%   iv. energy (E-E0) as a function of time
%   v. horizontal displacement of center of mass as a function of time

close all;  % close any open figures

% Initialize known values
m = 0.5;        % mass of each particle, kg
L = 5;          % rod length, m
g = 9.81;       % gravitation constant, m/s^2
theta_0 = 75;   % initial rod angle, deg
thetadot_0 = 0; % initial rod rotation rate, deg/s
x_0 = 1;        % initial block position, m
xdot_0 = 0;         % initial block velocity, m/s
t0 = 0;         % initial time, s
tf = 100;       % final time, s

% Initialize integration properties
dt = 0.01;      % time step, s
tnext = t0+dt;  % next time, s
options = odeset('RelTol', 1E-12, 'AbsTol', 1E-12, 'InitialStep', dt,...
    'MaxStep', dt);
c = 1;          % counter

% Initialize vectors for state
limit = tf/dt;
state = zeros(limit, 4);    % matrix to hold state variables
y0 = [x_0 deg2rad(theta_0) xdot_0 deg2rad(thetadot_0)].'; 
    % initial state vector
state(1,:) = y0;

% Initialize vectors for plots
x = zeros(1, limit+1);      % block position, m
x(1) = x_0;
x_dot = zeros(1, limit+1);  % block velocity, m/s
x_dot(1) = xdot_0;
theta = zeros(1, limit+1);  % rod angle, deg
theta(1) = theta_0;
theta_dot = zeros(1, limit+1);  % rod rotation rate, deg/s
theta_dot(1) = thetadot_0;
t = t0:dt:tf;     % time, s

% Numerically integrate for system
tic;    % start timer to check for excessively long loop
cycles = 0; % number of cycles so far
while tnext <= tf && cycles < 3
    c=c+1;  % increment counter
    [tn, yn] = ode45('block_rod_system', [t0 tnext], y0, options);
        % integrate
    % set up next
    t0 = tn(end);
    tnext = tnext + dt;
    y0 = yn(end, 1:4);
    state(c,:) = y0;    % save state at this time
    x(c) = y0(1);
    theta(c) = rad2deg(y0(2));
    x_dot(c) = y0(3);
    theta_dot(c) = rad2deg(y0(4));
    if (toc > 60*5)     % throw error if in loop longer than 5 minutes
        error('We''ve been looping for 5 minutes. Speed it up!');
    end
    if (abs(theta(c) - theta_0) < 0.00125) 
        % count number of cycles (theta(c) may not EXACTLY equal theta_0)
        cycles = cycles + 1;
    end        
end

% Print time for 3 cycles
fprintf('\nProblem 1a\ntime for 3 cycles\n');
disp(tnext-dt);

% Calculate energy and center of mass position
V = m*g*L*(cosd(75)-cosd(theta));
Ta = 1/2*m*x_dot.^2;
Tb = 1/2*m*((x_dot + L*deg2rad(theta_dot).*cosd(theta)).^2 ...
    + (L*deg2rad(theta_dot).*sind(theta)).^2);
E = V + Ta + Tb;
x_cm0 = x_0 + L/2*sind(theta_0);
x_cm = x + L/2*sind(theta) - x_cm0;

% Create and format plots
figure(1); plot(t(1:c), theta(1:c)); grid on;
title('Problem 1b, Plot i');xlabel('t');ylabel('theta(t)');
figure(2);plot(t(1:c), x(1:c));grid on;
title('Problem 1b, Plot ii');xlabel('t');ylabel('x(t)');
figure(3);plot(x(1:c), theta(1:c));grid on;
title('Problem 1b, Plot iii');xlabel('x');ylabel('theta(x)');
figure(4);plot(t(1:c), E(1:c));grid on;
title('Problem 1b, Plot iv');xlabel('t');ylabel('E - E_0');
figure(5);plot(t(1:c), x_cm(1:c));grid on;
title('Problem 1b, Plot v');xlabel('t');ylabel('x_{cm}(t)');


%%% Problem 2 %%%
% Part c
% Transform the given inertia matrix into the dashed system
I = [3 -2 1;-2 3 1;1 1 3];  % inertia matrix in original system
l = [0.766 0 0.643;0 0 0;-0.643 0 0.766];   % rotation matrix from
                                            % original system to
                                            % dashed system
Idash = l*I*l.';    % inertia matrix in dashed system
fprintf('\nProblem 2\ninertia matrix in dashed system\n');
disp(Idash);