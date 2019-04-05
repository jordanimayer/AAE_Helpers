%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 05
%
% Numerically simulate rotational motion of axisymmetric rigid body S
% moving in inertial frame N.
%%%%%

%% Preliminary setup

close all; clear all; format compact; rehash toolbox;

%% Problem 1b
% calculate results from analytical solutions for angular velocities

w0 = [2 1 -1];  % initial angular velocity vector, rad/s, expressed in
                % s_hat frame
                
% calculate results
t = [0:0.1:30].';  % span of time to simulate over, s
w1_an = ones(size(t,1),1)*w0(1);  % s_hat1 component, rad/s (constant)
w2_an = -sin(0.7*t) + 0.8429*cos(0.7*t) + 0.1571;
  % s_hat2 component, rad/s
w3_an = -cos(0.7*t) - 0.8429*sin(0.7*t);  % s_hat3 component, rad/s

% create plots
figure();
plot(t, w1_an, '-b'); hold on; plot(t, w2_an, '-r'); plot(t, w3_an, '-g');
title('Jordan Mayer - AAE 590/440, HW 05, Problem 1b: Analytical Results');
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); grid on;
ylim([-2.5 2.5]);
legend('s-hat_1 component', 's-hat_2 component', 's-hat_3 component');

t_an = t;  % for plotting along with numerical results, s
                
%% Problem 1c,d
% numerically simulate dynamic and kinematic differential equations
% simultaneously

% set up initial conditions
q_0 = [0 0 0 1];  % initial Euler parameters, [vector scalar]
                  % vector components expressed in s_hat frame
w_0 = [2 1 -1];  % initial angular velocities, rad/s
                 % measure numbers expressed in s_hat frame
y_0 = [q_0 w_0].';  % initial state vector

% numerically integrate
t_0 = 0;  % initial time, s
t_f = 30;  % final time, s
t_span = [t_0 t_f];  % span of time to integrate over
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
  % error tolerance is set to 10^(-7)
[t, y] = ode45('hw05_numint_qw', t_span, y_0, options);

% extract results
q1 = y(:,1);  % s_hat1 Euler parameter
q2 = y(:,2);  % s_hat2 Euler parameter
q3 = y(:,3);  % s_hat3 Euler parameter
q4 = y(:,4);  % scalar Euler parameter
K = q1.^2 + q2.^2 + q3.^2 + q4.^2;
  % Euler parameters constraint; should be 1
K_0 = K(1);  % initial value for K
errK = K - K_0;  % error for K
w1_num = y(:,5);  % s_hat1 angular velocity, rad/s
w2_num = y(:,6);  % s_hat2 angular velocity, rad/s
w3_num = y(:,7);  % s_hat3 angular velocity, rad/s
t_num = t;  % for plotting along with analytical results, s

% create plots
figure();
plot(t, errK, '-b'); grid on;
title('Jordan Mayer - AAE 590/440, HW 05, Problem 1c: Numerical Results');
xlabel('Time (s)'); ylabel('K - K_0');
figure();
plot(t, w1_num, '-b'); hold on; plot(t, w2_num, '-r');
plot(t, w3_num, '-g');
title('Jordan Mayer - AAE 590/440, HW 05, Problem 1d: Numerical Results');
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); grid on;
ylim([-2.5 2.5]);
legend('s-hat_1 component', 's-hat_2 component', 's-hat_3 component');

% plot together with analytical
figure();
plot(t_an, w1_an, '-b'); hold on; plot(t_an, w2_an, '-r'); plot(t_an, w3_an, '-g');
plot(t_num, w1_num, '--b'); plot(t_num, w2_num, '--r'); plot(t_num, w3_num, '--g');
title('Jordan Mayer - AAE 590/440, HW 05, Problem 1b: Both Angle Results');
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); grid on;
ylim([-2.5 2.5]);
legend('Analytical w_1', 'Analytical w_2', 'Analytica w_3', ...
       'Numerical w_1', 'Numerical w_2', 'Numerical w_3');

% calculate final value differences
format long;  % necessary to see difference
wf_an = [w1_an(end) w2_an(end) w3_an(end)]
wf_num = [w1_num(end) w2_num(end) w3_num(end)] 
wf_diff = wf_an - wf_num % rad/s
format short;

%% Problems 2b and 3

% set initial conditions (vectors in s_hat frame)
q_0 = [0 0 0 1];
w_0 = [2 1 -1];  % rad/s
C_0 = eye(3);  % initial direction cosine matrix from N to S
y_0 = [q_0 w_0 C_0(1,:) C_0(2,:) C_0(3,:)].';  % initial state vector

% numerically integrate
t_0 = 0;  % initial time, s
t_f = 30;  % final time, s
t_span = [t_0 t_f];  % span of time to integrate over
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
  % error tolerance is set to 10^(-7
[t, y] = ode45('hw05_numint_all', t_span, y_0, options);
  % default relative error tolerance for integration is 1e-3 (0.001)
  
% print output headers
fprintf('\t t (s) \t w_1 (rad/s) \t w_2 (rad/s) \t w_3 (rad/s) \t C_11 \t C_21 \t C_31 \t beta (deg) \t sigma (deg) \t (K-K_0)*10^9 \n');
fprintf('\t-------\t-------------\t-------------\t-------------\t------\t------\t------\t-------------\t------------\t---------\n');
  
% extract states (vectors in s_hat frame)
q1 = y(:,1);
q2 = y(:,2);
q3 = y(:,3);
q4 = y(:,4);
w1 = y(:,5);  % rad/s
w2 = y(:,6);  % rad/s
w3 = y(:,7);  % rad/s
C11 = y(:,8);
C21 = y(:,11);
C31 = y(:,14);
K = q1.^2 + q2.^2 + q3.^2 + q4.^2;  % constraint constant, should be 1
K0 = K(1);
errK = K - K0;

n = size(y,1);  % number of output stats
beta = zeros(n,1);  % precession angle, deg
sigma = zeros(n,1);  % nutation angle, deg
beta(1) = 0; sigma(1) = 0;  % initial conditions
                            % (avoid divide by zero in for loop)

k = 1;
fprintf('\t %.2f \t %.2f \t\t\t %.2f \t\t\t %.2f \t\t\t %.2f \t %.2f \t %.2f \t %.2f   \t\t %.2f \t\t\t %.2f \n',...
                    t(k), w1(k), w2(k), w3(k), C11(k), C21(k), C31(k), beta(k), sigma(k), errK(k)*10^9);
prev_t = 0;  % for printing outputs (see below)
for k=2:n
    
    % get nutation angle
    sigma_k = acosd(C11(k));  % deg
    % adjust to convention (0deg < sigma < 180deg)
    if sigma_k < -180
        sigma_k = 360 + sigma_k;
    elseif sigma_k < 0
        sigma_k = -sigma_k;
    elseif sigma_k >= 180
        sigma_k = -(sigma_k - 360);
    end
    if (sigma_k < 0 || sigma_k > 180)
        error('sigma_k is not b/w 0 and 180!');
    end
    
    % get precession angle
    beta_k = zeros(2,2);  % top row: 2 options from inverse sine
                          % bottom row: 2 options from inverse cosine
    beta_k(1,1) = asind(C31(k)/sind(sigma_k));  % deg
    beta_k(1,2) = 180 - beta_k(1,1);  % deg
    beta_k(2,1) = acosd(C21(k)/sind(sigma_k));  % deg
    beta_k(2,2) = -beta_k(2,1);  % deg
    % ensure between -180 and 180 (for comparison and plotting)
    for e=1:4
        while (beta_k(e) < -180)
            beta_k(e) = beta_k(e) + 360;
        end
        while (beta_k(e) > 180)
            beta_k(e) = beta_k(e) - 360;
        end
    end
    % choose value with least difference between inverse sine result and
    % inverse cosine result
    diff_11 = abs(beta_k(1,1) - beta_k(2,1));
    diff_12 = abs(beta_k(1,1) - beta_k(2,2));
    diff_21 = abs(beta_k(1,2) - beta_k(2,1));
    diff_22 = abs(beta_k(1,2) - beta_k(2,2));
    diff_min = min([diff_11 diff_12 diff_21 diff_22]);
    if diff_11 == diff_min || diff_12 == diff_min
        beta_k = beta_k(1,1);  % deg
    elseif diff_21 == diff_min || diff_22 == diff_min
        beta_k = beta_k(1,2);  % deg
    else
        disp(diff_11); disp(diff_12); disp(diff_21); disp(diff_22); disp(beta_k);
        error('no matching values for beta_k!');
    end
    
    sigma(k) = sigma_k;
    beta(k) = beta_k;
    
    % print output and ~2 second intervals
    if mod(t(k),2) < 0.02
        % t is close to a multiple of 2
        if (t(k) - prev_t >= 1.99)
            % not too close to previous time
            prev_t = t(k);
            fprintf('\t %.2f \t %.2f \t\t\t %.2f \t\t\t %.2f \t\t\t %.2f \t %.2f \t %.2f \t %.2f   \t\t %.2f \t\t\t %.2f \n',...
                    t(k), w1(k), w2(k), w3(k), C11(k), C21(k), C31(k), beta(k), sigma(k), errK(k)*10^9);
        end
    end
end

% create plots for Problem 3
figure();
plot(t, beta, '-b'); hold on; plot(t, sigma, '-r');
title('Jordan Mayer - AAE 590/440, HW 05, Problem 3a: Numerical Results');
xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;
legend('Precession', 'Nutation');
figure();
plot(C21, C31, '-b');
title('Jordan Mayer - AAE 590/440, HW 05, Problem 3b: Numerical Results');
xlabel('C_{21}'); ylabel('C_{31}'); grid on;
axis square; axis([-0.4 1 -1 0.4]);

%% Problem 4
% will change value of K2 in hw05_numint_all and then run this

% set initial conditions (vectors in s_hat frame)
q_0 = [0 0 0 1];
w_0 = [2 1 -1];  % rad/s
C_0 = eye(3);  % initial direction cosine matrix from N to S
y_0 = [q_0 w_0 C_0(1,:) C_0(2,:) C_0(3,:)].';  % initial state vector

% numerically integrate
t_0 = 0;  % initial time, s
t_f = 30;  % final time, s
t_span = [t_0 t_f];  % span of time to integrate over
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
  % error tolerance is set to 10^(-7
[t, y] = ode45('hw05_numint_all', t_span, y_0, options);
  % default relative error tolerance for integration is 1e-3 (0.001)

% extract states
C11 = y(:,8);
C21 = y(:,11);
C31 = y(:,14);

n = size(y,1);  % number of output stats
sigma = zeros(n,1);  % nutation angle, deg
sigma(1) = 0;  % initial condition

for k=2:n
    
    % get nutation angle
    sigma_k = acosd(C11(k));  % deg
    % adjust to convention (0deg < sigma < 180deg)
    if sigma_k < -180
        sigma_k = 360 + sigma_k;
    elseif sigma_k < 0
        sigma_k = -sigma_k;
    elseif sigma_k >= 180
        sigma_k = -(sigma_k - 360);
    end
    if (sigma_k < 0 || sigma_k > 180)
        error('sigma_k is not b/w 0 and 180!');
    end
    
    sigma(k) = sigma_k;
end

% create plots
figure();
plot(t, sigma, '-b');
title('Jordan Mayer - AAE 590/440, HW 05, Problem 3a: Numerical Results, T = 0 N-m');
xlabel('Time (s)'); ylabel('Nutation Angle (deg)'); grid on;
figure();
plot(C21, C31, '-b');
title('Jordan Mayer - AAE 590/440, HW 05, Problem 3b: Numerical Results, T = 0 N-m');
xlabel('C_{21}'); ylabel('C_{31}'); grid on;
axis square; axis([-0.4 1 -1 0.4]);