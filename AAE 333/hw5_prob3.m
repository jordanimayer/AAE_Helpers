% AAE 333: Fluide Mechanics - HW 5, Problem 3
%
% This program solves a differential equation to find the velocity of a
% rocket

t0 = 0.01;     % Start time
tf = 1.56;     % End time (after this the rocket is not accelerating and
               % mass is constant, so the diffEQ no longer holds.

tspan = t0:0.01:tf; % Span of time

v0 = 0;     % Initial velocity

% The diffEQ to be solved
dvdt = @(t, v) 1957/(10-5.76*t) - 0.00146*(v^2)/(10-5.76*t) - 9.81;

[t, v] = ode45(dvdt, tspan, v0);    % Solve diffEQ

% Plot results

% Without drag
v_t = -340*log(abs(t-1.74)) - 0.81*t + 188;
plot(t, v_t, '^k');
xlabel('Time, s');
ylabel('Velocity, m/s');
title('Small Rocket Acceleration Prediction');
grid on;
hold on;

% With drag
plot(t, v, '-k');
xlabel('Time, s');
ylabel('Velocity, m/s');
title('Small Rocket Acceleration Prediction');
grid on;
legend('Without Drag', 'With Drag');

% Display results
fprintf('Maximum velocity without drag at time t= 1.56 s: 756 m/s\n');
fprintf('Maxium velocity with drag at time %.2f: %.0f m/s\n', 1.56, v(find(t==1.56, 1)));