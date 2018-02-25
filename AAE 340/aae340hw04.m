%%%%%
% AAE 340 - HW 04
% Jordan Mayer
%%%%%

%%% Problem 2d %%%
% Numerically integrate EOMs for satellite in Keplerian orbit over
% one orbital period. Initial conditions: geostationary ring,
% completing one revolution around center of Earth in 24h,
% at distance of 42,000 km from Earth center.
% Plot result in 3D
%
% NOTE: assume r_dot0 = 0;
%       Keplerian --> treat z as constant = 0;

% Known values
mu = 0.39860*10^6*(10^3)^3;     % mu (GM) for Earth (m^3/s^2)
                                % source: NASA GSFC website
T = 24*60*60;                   % orbital period (s)
t0 = 0;                         % initial time (s)
r0 = 42000*10^3;                % initial radius (m)
r_dot0 = 0;                     % initial radial velocity (m/s)
theta0 = 0;                     % initial angle (rad)
theta_dot0 = sqrt(mu/r0^3);     % initial angular velocity (rad/s)
                                % source: AAE 251 notes

% Initialize integration properties
dt = 1.0;                       % time step (s)
tnext = t0+dt;                  % next time (s)
options = odeset('RelTol',1E-12,'AbsTol',1E-12,'InitialStep',dt,...
    'MaxStep',dt);
c=1;    % counter

% Initialize vectors for state
limit = T/dt;
state = zeros(limit, 4);    % matrix to hold state variables
y0 = [r0 r_dot0 theta0 theta_dot0]';    % initial state vector
state(1,:) = y0;

% Initialize x/y/z for 3D plot
x = zeros(1,limit+1);
y = zeros(1,limit+1);
z = zeros(1,limit+1);
x(1) = r0*cos(theta0);
y(1) = r0*sin(theta0);

% Numerically integrate for orbit
tic;    % start timer to check for infinite loop
while tnext <= T
    c=c+1;  % increment counter
    [tn,yn] = ode45('orbitpde',[t0 tnext],y0,options);  % integrate
    % set up next time
    t0=tn(end);
    tnext=tnext+dt;
    y0=yn(end,1:4);
    state(c,:)=y0;  % save state at this time
    % calculate x/y for 3D plot
    r=y0(1);        % radius (m)
    theta=y0(3);    % theta (rad)
    x(c) = r*cos(theta);
    y(c) = r*sin(theta);
    if (toc > 60*5)     % throw error if in loop longer than 5 minutes
        error('Wait, stop. There might be an infinite loop here.');
    end
end

% Create and format 3D plot of orbit
plot3(x,y,z,'-r');
xlabel('x (m)');ylabel('y (m)');zlabel('z (m)');
title('Keplerian Orbit of Satellite');grid on;