%%%%%
% AAE 340 - HW 07
% Jordan Mayer
%%%%%

%%% Problem 2c-2f %%%
% Numerically integrate EOMs for binary star system over 1 billion seconds.
% State space EOMs derived in parts a and b.
% Initial conditions given.
%
% Problem 2c
% Plot positions of both stars relative to origin of inertial frame.
%
% Problem 2d
% Plot position of center of mass.
%
% Problem 2e
% Plot velocity of center of mass and explain
% (explanation attached with plot).
%
% Problem 2f
% Plot positions of both stars relative to center of mass.

close all;      % close any open figures

% Initialize known values
G = 6.673*10^(-11);         % universal graviational constant, m^3/(kg*s)
m_sun = 1.987*10^30;        % mass of sun, kg
r1_0 = [0, 0];              % initial position of star 1, km
r2_0 = [0, 1.000*10^9];     % initial position of star 2, km
rdot1_0 = [2.000, 0];       % initial velocity of star 1, km/s
rdot2_0 = [12.000, 4.000];  % initial velocity of star 2, km/s
t0 = 0;                     % initial time, s
tf = 1.000*10^9;            % final time, s
    
% Calculate m1, m2
m1 = m_sun;                 % mass of star 1, kg
m2 = m1/2;                  % mass of star 2, kg

% Convert r and rdot to x, y, xdot, and ydot
x1_0 = r1_0(1)*10^3;        % initial x-position of star 1, m
y1_0 = r1_0(2)*10^3;        % initial y-position of star 1, m
x2_0 = r2_0(1)*10^3;        % initial x-position of star 2, m
y2_0 = r2_0(2)*10^3;        % initial y-position of star 2, m
xdot1_0 = rdot1_0(1)*10^3;  % initial x-velocity of star 1, m/s
ydot1_0 = rdot1_0(2)*10^3;  % initial y-velocity of star 1, m/s
xdot2_0 = rdot2_0(1)*10^3;  % initial x-velocity of star 2, m/s
ydot2_0 = rdot2_0(2)*10^3;  % initial y-velocity of star 2, m/s

% Initialize integration properties
dt = 10^4;                  % time step, s
                            % (use smallest time step possible while 
                            % staying below 5 minutes of computation time)
tnext = t0+dt;              % next time, s
options = odeset('RelTol', 1E-12, 'AbsTol', 1E-12, 'InitialStep', dt,...
    'MaxStep', dt);         
c = 1;                      % counter

% Initialize vectors for state
limit = tf/dt;
state = zeros(limit, 8);    % matrix to hold state variables
y0 = [x1_0 y1_0 x2_0 y2_0 xdot1_0 ydot1_0 xdot2_0 ydot2_0]';
                            % initial state vector
state(1,:) = y0;

% Initialize vectors for plots
x1 = zeros(1, limit+1);     % x-position of star 1, m
y1 = zeros(1, limit+1);     % y-position of star 1, m
x2 = zeros(1, limit+1);     % x-position of star 2, m
y2 = zeros(1, limit+1);     % y-position of star 2, m
xdot1 = zeros(1, limit+1);  % x-velocity of star 1, m
ydot1 = zeros(1, limit+1);  % y-velocity of star 1, m
xdot2 = zeros(1, limit+1);  % x-velocity of star 2, m
ydot2 = zeros(1, limit+1);  % y-velocity of star 2, m
x1(1) = x1_0;
y1(1) = y1_0;
x2(1) = x2_0;
y2(1) = y2_0;
xdot1(1) = xdot1_0;
ydot1(1) = ydot1_0;
xdot2(1) = xdot2_0;
ydot2(1) = ydot2_0;

% Numerically integrate for orbit
tic;    % start timer to check for infinite loop
while tnext <= tf
    c = c+1;  % increment counter
    [tn, yn] = ode45('binarystar', [t0 tnext], y0, options);  % integrate
    % set up next time
    t0 = tn(end);
    tnext = tnext + dt;
    y0 = yn(end, 1:8);
    state(c, :) = y0;   % save state at this time
    x1(c) = y0(1);
    y1(c) = y0(2);
    x2(c) = y0(3);
    y2(c) = y0(4);
    if (toc > 60*5)     % throw error if in loop longer than 5 minutes
        error('We''ve been looping for 5 minutes. Speed it up!');
    end
end

% Calculate values related to center of mass

xc = 1/(m1+m2) * (m1*x1 + m2*x2);   % x-position of center of mass, m
yc = 1/(m1+m2) * (m1*y1 + m2*y2);   % y-position of center of mass, m
xdotc = 1/(m1+m2) * (m1*xdot1 + m2*xdot2);  
                                  	% x-velocity of center of mass, m
ydotc = 1/(m1+m2) * (m1*ydot1 + m2*ydot2);  
                                    % y-velocity of center of mass, m
xcp1 = x1 - xc;   % x-position of star 1 relative to COM, m
ycp1 = y1 - yc;   % y-position of star 1 relative to COM, m
xcp2 = x2 - xc;   % x-position of star 2 relative to COM, m
ycp2 = y2 - yc;   % y-position of star 2 relative to COM, m

% Create and format plots

figure(1);              % plot positions of both stars relative to
                        % origin of inertial frame
plot(x1/1000, y1/1000, '-b'); hold on; 
plot(x2/1000, y2/1000, '-r'); grid on;
axis([0 6*10^9 0 6*10^9]); pbaspect([1 1 1]);
    % scale axes for easy comparison and interpretation
title('Inertial Motion of Binary Star System over 1 Billion Seconds');
xlabel('x-position, km'); ylabel('y-position, km');
legend('Star 1', 'Star 2');
hold off;

figure(2);              % plot position of center of mass
plot(xc/1000, yc/1000, '-k'); grid on;
axis([0 6*10^9 0 6*10^9]); pbaspect([1 1 1]);
title('Motion of Binary Star System COM over 1 Billion Seconds');
xlabel('x-position, km'); ylabel('y-position, km');

figure(3);              % plot velocity of center of mass
plot(xdotc/1000, ydotc/1000, 'ok'); grid on;
axis([0 6 0 6]); pbaspect([1 1 1]);
title('Velocity of Binary Star System COM over 1 Billion Seconds');
xlabel('x-velocity, km/s'); ylabel('y-velocity, km/s');

figure (4);             % plot positions of both stars relative to
                        % center of mass
plot(xcp1/1000, ycp1/1000, '-b'); hold on;
plot(xcp2/1000, ycp2/1000, '-r'); grid on;
axis([-4*10^8 8*10^8 -4*10^8 8*10^8]); pbaspect([1 1 1]);
title('Relative Motion of Binary Star System over 1 Billion Seconds');
xlabel('x-position (from COM), km'); ylabel('y-position (from COM), km');
legend('Star 1', 'Star 2');