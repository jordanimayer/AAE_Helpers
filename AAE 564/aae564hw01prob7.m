%%%%%
% AAE 564
% HW 01
% Jordan Mayer
%%%%%

%%% Problem 7 %%%
% Simulate double-pendulum system for various parameters and initial
% conditions.
% EOMs in class notes. State space description in handwritten homework.

close all;  % close any open figures

% Initialize initial conditions [IC1 IC2 IC3 IC4 IC5 IC6 IC7]
y_ICs = zeros(1,7);
theta1_ICs = [-10 10 -90 -90.01 100 100.01 179.99];  % deg
theta2_ICs = [10 10 90 90 100 100 0];  % deg
y_dot_ICs = zeros(1,7);
theta1_dot_ICs = zeros(1,7);
theta2_dot_ICs = zeros(1,7);

P = 4;  % CHANGE THIS to change parameter set

% Simulate for the following combinations:
% P1: IC1, IC2, IC3, IC7
% P4: IC1, IC2, IC3, IC4
tic;  % start timer to check for excessively long loop
IC = zeros(1,4);
if (P == 1)
    IC = [1 2 3 7];
elseif (P == 4)
    IC = [1 2 3 4];
end
    
for k = IC
    % Initialize integration properties
    t0 = 0;  % initial time
    tf = 50;  % final time
    dt = 0.01;  % time step
    tnext = t0 + dt;  % next time
    options = odeset('RelTol', 1E-12, 'AbsTol', 1E-12, 'InitialStep', dt,...
        'MaxStep', dt);
    c = 1;  % counter

    % Initialize state vectors
    limit = tf / dt;
    state = zeros(limit, 6);  % matrix to hold state variables

    % Initialize vectors for plots
    y = zeros(1, limit+1);
    y_dot = zeros(1, limit+1);
    theta1 = zeros(1, limit+1);
    theta1_dot = zeros(1, limit+1);
    theta2 = zeros(1, limit+1);
    theta2_dot = zeros(1, limit+1);
    t = t0:dt:tf;

    % Set up initial conditions
    y_0 = y_ICs(k);
    y(1) = y_0;
    theta1_0 = theta1_ICs(k);
    theta1(1) = theta1_0;
    theta2_0 = theta2_ICs(k);
    theta2(1) = theta2_0;
    y_dot_0 = y_dot_ICs(k);
    y_dot(1) = y_dot_0;
    theta1_dot_0 = theta1_dot_ICs(k);
    theta1_dot(1) = theta1_dot_0;
    theta2_dot_0 = theta2_dot_ICs(k);
    theta2_dot(1) = theta2_dot_0;
    x0 = [y_0, y_dot_0, theta1_0, theta1_dot_0, theta2_0, theta2_dot_0].';
      % initial state vector
    
    while tnext <= tf
        c = c+1;  % increment counter
        
        [tn, xn] = ode45('double_pendulum', [t0 tnext], x0, options);
          % integrate state vector
        
        % set up next
        t0 = tn(end);
        tnext = tnext + dt;
        x0 = xn(end, 1:6);
        state(c,:) = x0;  % save state at this time
        y(c) = x0(1);
        y_dot(c) = x0(2);
        theta1(c) = x0(3);
        theta1_dot(c) = x0(4);
        theta2(c) = x0(5);
        theta2_dot(c) = x0(6);
        
        % check for loop longer than 5 minutes
        if (toc > 60*5)
            error('We''ve been looping for 5 minutes. Speed it up!');
        end
    end
    
    % Create and format plots
    ttl = ['Jordan Mayer: Problem 7 Results'...
           ' P' num2str(P) ' IC' num2str(k)];
    figure('pos', [50 50 900 900]);
    subplot(3,2,1); plot(t, y); grid on; title(ttl); xlabel('t');
    ylabel('y(t)');
    subplot(3,2,2); plot(t, y_dot); grid on; title(ttl); xlabel('t');
    ylabel('ydot(t)');
    subplot(3,2,3); plot(t, theta1); grid on; title(ttl); xlabel('t');
    ylabel('theta1(t) (deg)');
    subplot(3,2,4); plot(t, theta1_dot); grid on; title(ttl); xlabel('t');
    ylabel('theta1dot(t) (deg/time)');
    subplot(3,2,5); plot(t, theta2); grid on; title(ttl); xlabel('t');
    ylabel('theta2(t) (deg)');
    subplot(3,2,6); plot(t, theta2_dot); grid on; title(ttl); xlabel('t');
    ylabel('theta2dot(t) (deg/time)');
end