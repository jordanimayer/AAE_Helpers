%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 09
%
% hw09_sim_and_plot:
%   Simulate body motion using numerical integration and create plots of
%   nutation angle and error.
%
% Inputs:
%   v_f: number of orbital revolutions to simulate
%   sigma_0: initial nutation angle, deg
%   r: initial angular velocity about orbit normal axis, divided by orbital
%      rate
%   next_fig: next figure number
%   plot_err: 1 to plot error, 0 otherwise
%
% Outputs:
%   None
%%%%%

function [] = hw09_sim_and_plot(v_f, sigma_0, r, next_fig, plot_err)
    % set up constants
    v_0 = 0;
    v_span = [v_0, v_f];
    options = odeset('RelTol', 1.0e-7, 'AbsTol', 1.0e-7);
      % set error tolerance to 10^(-7)
    sigma_max_diff = 1;
    
    % set up initial state
    w_0 = [r 0 0];  % nondimensionalized angular velocities of S in N,
                    % expressed in c_hat frame
    q_0 = [0 0 sind(sigma_0/2) cosd(sigma_0/2)];
      % quaternion from A to C, vector in c_hat basis
    K_0 = sum(q_0.^2);  % quaternion constraint, should be 1.0 always
    y_0 = [w_0 q_0].';  % initial state vector
    [v, y] = ode45('hw09_numint_all', v_span, y_0, options);
    
    % extract states
    w1 = y(:,1);
    w2 = y(:,2);
    w3 = y(:,3);
    q1 = y(:,4);
    q2 = y(:,5);
    q3 = y(:,6);
    q4 = y(:,7);
    
    % calculate nutation angles
    sigma = zeros(size(y,1),1);
    sigma(1) = sigma_0;
    for k = 2:size(y,1)
        sigma(k) = theta_from_quaternion([q1(k) q2(k) q3(k) q4(k)], ...
                                         sigma_max_diff);
    end
    
    % calculate K values
    K = (q1.^2 + q2.^2 + q3.^2 + q4.^2).^(1/2);
    deltaK = K - K_0;
    
    % plot nutation angle
    figure(next_fig);
    plot(v, sigma); hold on;
    if plot_err == 1
        % plot error, as measured by quaternion constraint
        figure(next_fig+1);
        plot(v, deltaK); hold on;
    end
end