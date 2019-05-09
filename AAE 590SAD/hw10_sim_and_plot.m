%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 10
%
% hw10_sim_and_plot:
%   Simulate body motion using numerical integration and create plots of
%   three angles: nutation (sigma), precession relative to orbit (eta), and
%   inertial precession (lambda)
%
% Inputs:
%   v_span: span of orbital revolutions to simulate, [start, stop]
%   sigma_0: initial nutation angle, deg
%   r: initial angular velocity about orbit normal axis, divided by orbital
%      rate
%   next_fig: next figure number
%
% Outputs:
%   None
%%%%%

function [] = hw10_sim_and_plot(v_span, sigma_0, r, next_fig)
    % set constants
    J = 150;  % principal moment of inertia about body axis of symmetry,
              % kg*m^2
    I = 500;  % principal moment of inertia about transverse body axes,
              % kg*m^2
    
    % set up integration
    options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10);
      % set error tolerance to 10^(-10)
    angle_max_diff = 1;  % maximum difference for two angle results to be
                         % considered equal (used for quadrant checking)
                         
    % set up initial state and integrate
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
    
    % calculate angles
    sigma = zeros(size(y,1),1); eta = sigma; lambda = sigma;
    lambda_part = sigma;
    delta_lambda = 3*pi/r * (I/J - 1) * cosd(sigma_0);
      % rate of change of nominal lambda (rad/rev)
    factor360 = 0;
    for k = 1:size(y,1)
        C = dircos_from_quaternion([q1(k), q2(k) q3(k) q4(k)]);
          % direction cosine matrix from A to C
        C11 = C(1,1); C21 = C(2,1); C31 = C(3,1);
        sigma(k) = acosd(C11);
        eta(k) = theta_from_sin_cos(C31/sind(sigma(k)), ...
                                    C21/sind(sigma(k)), angle_max_diff);
        if (k > 1) && ((eta(k-1) - eta(k)) > 180)
            % eta has "flipped" from 360 to 0 degrees
            factor360 = factor360 + 1;
        elseif (k > 1) && ((eta(k-1) - eta(k)) < -180)
            % eta has "flipped" from 0 to 360 degrees
            factor360 = factor360 - 1;
        end
        gamma = rad2deg(2*pi*v(k));  % inertial precession of orbit, deg
        lambda(k) = gamma + eta(k) + 360*factor360;
        lambda_part(k) = rad2deg(delta_lambda*v(k));
          % lambda estimate using delta_lambda
    end
    
    % plot angles
    figure(next_fig);
    plot(v, sigma); hold on;
    figure(next_fig+1);
    plot(v, lambda); hold on;
      % comment-out if plotting both lambda and delta_lambda estimation
%     figure(next_fig+1);
%     plot(v, eta); hold on;
%     if r ~= 200
%         c = '';  % color of lambda plots for this r value
%         switch abs(r)
%             case 1.5
%                 c = 'r';
%             case 8
%                 c = 'm';
%             case 22
%                 c = 'g';
%             case 120
%                 c = 'c';
%             otherwise
%                 error('No lambda color for this r value!');
%         end
%         figure(next_fig+2);
%         plot(v, lambda, ['-', c]); hold on;
%         plot(v_span, [lambda(1), lambda(end)], ['--', c]);
%     end  % comment-out if plotting only lambda
end