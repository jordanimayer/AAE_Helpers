%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 11
%
% hw11p4_sim_and_plot:
%   Simulate body motion using numerical integration and create plots of
%   three angles: nutation (sigma), precession relative to orbit (eta), and
%   inertial precession (lambda).
%   Includes rotor (for use in Problem 4).
%
% Inputs:
%   v_span: span of orbital revolutions to simulate, [start, stop]
%   region: which region of the K3-K2 plot to simulate for
%   w_0: initial nondimensionalized angular velocity vector (B in N)
%   q_0: initial Euler parameters (A to B)
%
% Outputs:
%   None
%%%%%

function [] = hw11p4_sim_and_plot(v_span, region, w_0, q_0)
    global I_reg1 I_reg4 I_reg6 I n_rotor
    % I_reg1, I_reg4, I_reg6: gyrostat principal moments of inertia in
    %                         terms of body axes, kg-met^2, for each region
    % I: for specific region we're analyzing, for use in integration
    % n_rotor: nondimensionalized rotor spin rate
    
    
    % set plot formatting variables  
    ttl_start = ['Jordan Mayer, AAE 590/440, PS11, Problem 4, ', ...
                 'Numerical Results, '];
    ttl = [ttl_start, 'n = ', num2str(n_rotor), ...
           ', Region ', num2str(region)];
    lgd = {'particular solution', ...
           'w_{20} = w_{30} = 0.08', 'w_{20} = w_{30} = 0.16'};
    xlab = 'Number of Revolutions';
    fpos = [9 49 944 918];
    
    % set shape variables, figure numbers
        % gyrostat principal moments of inertia in terms of body axes,
        % kg-met^2
    if region == 1
        % orientation ii
        I = I_reg1;
        fnum = 1;
    elseif region == 4
        % orientation i
        I = I_reg4;
        fnum = 2;
    elseif region == 6
        % orientation iii
        I = I_reg6;
        fnum = 3;
    end
    
    if w_0(1) == 1.0 && w_0(2) == 0 && w_0(3) == 0
        % particular solution - no need for integration!
        v = v_span;
        sigma = [0, 0];
        lambda = [0, rad2deg(2*pi*v_span(2))];
    else
        % set up integration
        options = odeset('RelTol', 1.0e-10, 'AbsTol', 1.0e-10);
          % set error tolerance to 10^(-10)
        angle_max_diff = 1;  % maximum difference for two angle results to be
                             % considered equal (used for quadrant checking)

        % set up initial state and integrate
        y_0 = [w_0 q_0].';  % initial state vector
        [v, y] = ode45('hw11p4_numint_all', v_span, y_0, options);

        % extract quaternions
        q1 = y(:,4);
        q2 = y(:,5);
        q3 = y(:,6);
        q4 = y(:,7);

        % calculate angles
        sigma = zeros(size(y,1),1); eta = sigma; lambda = sigma;
        cycles_sigma = 0; cycles_eta = 0;
          % indicates how many "cycles" this angle has undergone (how many
          % times to add 360 deg to asin/acos result)
        for k = 2:size(y,1)
            % start at second element, since all angles will be zero
            % initially
            C = dircos_from_quaternion([q1(k), q2(k) q3(k) q4(k)]);
              % direction cosine matrix from A to B
            C11 = C(1,1); C21 = C(2,1); C31 = C(3,1);

            sigma_k = acosd(C11);
            if (k > 1) && ((sigma(k-1) - sigma_k) > 180)
                % sigma has "flipped" from 360 to 0 degrees
                cycles_sigma = cycles_sigma + 1;
                %fprintf('cycle_sigma+ at k ='); disp(k);
            elseif (k > 1) && ((sigma(k-1) - sigma_k) < -180)
                % sigma has "flipped" from 0 to 360 degrees
                cycles_sigma = cycles_sigma - 1;
                %fprintf('cycle_sigma- at k ='); disp(k);
            end
            sigma(k) = sigma_k + 360*cycles_sigma;
              % nutation of body, deg

            eta_prev = eta(k-1) - 360*cycles_eta;
%             fprintf(['\n\nk = ', num2str(k)]);
%             fprintf(['\nC21 = ', num2str(C21)]);
%             fprintf(['\nC31 = ', num2str(C31)]);
            eta_k = theta_from_sin_cos(C31/sind(sigma_k), ...
                                       C21/sind(sigma_k), angle_max_diff);
            if (k > 1) && ((eta_prev - eta_k) > 180)
                % sigma has "flipped" from 360 to 0 degrees
                cycles_eta = cycles_eta + 1;
                %fprintf('cycle_eta+ at k ='); disp(k);
            elseif (k > 1) && ((eta_prev - eta_k) < -180)
                % sigma has "flipped" from 0 to 360 degrees
                cycles_eta = cycles_eta - 1;
                %fprintf('cycle_eta- at k ='); disp(k);
            end
            eta(k) = eta_k + 360*cycles_eta;
              % orbital precession of body, deg

            gamma_k = rad2deg(2*pi*v(k));
              % inertial precession of orbit, deg
            lambda(k) = gamma_k + eta(k);
              % inertial precession of body, deg
        end
    end
    
    fig = figure(fnum); set(fig, 'Position', fpos);
    
    % plot nutation angle
    subplot(2,1,1); plot(v, sigma); hold on;
    if w_0(2) == 0.16
        % last run; add formatting
        xlabel(xlab); ylabel('Nutation Angle (sigma), deg');
        title(ttl); grid on; legend(lgd);
    end
    
%     % plot orbital precession angle
%     subplot(3,1,2); plot(v, eta); hold on;
%     if w_0(2) == 0.16
%         % last run; add formatting
%         xlabel(xlab); ylabel('Orbital Precession Angle (sigma), deg'); grid on;
%     end
    
    % plot inertial precession angle
    subplot(2,1,2); plot(v, lambda); hold on;
    if w_0(2) == 0.16
        % last run; add formatting
        xlabel(xlab); ylabel('Inertial Precession Angle (sigma), deg');
        grid on;
    end
end