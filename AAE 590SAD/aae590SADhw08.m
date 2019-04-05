%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 08
%
% Perform numerical simulations and plot results.
% Perform basic arithmetic to obtain and verify other results.
%%%%%

%% Preliminary setup

close all; clear all; format compact; format short; rehash toolbox;

%% Problem 1
fprintf('\nProblem 1:\n');

fprintf('\nPart a:\n\n');

I_S_Sstar = [150 500 500];  % principal moments of inertia, kg*m^2

% calculate semi-diameters on corresponding axes
alpha = 1./sqrt(I_S_Sstar)
scale = 0.09;

% plot s1-s3 projection
x1 = linspace(-alpha(1),alpha(1),1000);
x3 = sqrt(1/I_S_Sstar(3) * (1 - I_S_Sstar(1)*x1.^2));
figure();
plot(x1, x3, '-b'); hold on; plot(x1, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 08, Problem 1a: Inertia Ellipsoid Projection (s1-s3)')
xlabel('s_1'); ylabel('s_3'); grid on; hold on;
axis square; axis([-1 1 -1 1] * scale);

% calculate angular velocity and momentum
N_omega_S = [2*cosd(45), 0, 2*sind(45)]  % rad/s, in s_hat frame
N_H_S_Sstar = [I_S_Sstar(1)*N_omega_S(1), 0, I_S_Sstar(3)*N_omega_S(3)]
  % kg*m^2/s, in s_hat frame
angle = atand(N_H_S_Sstar(3)/N_H_S_Sstar(1))
  % angle of N_H_S_Sstar with s_hat_1, deg
  
fprintf('\nPart b:\n\n');

% calculate precession rate, spin rate, and nutation angle
H = norm(N_H_S_Sstar);  % kg*m^2/s
I = I_S_Sstar(2);  % kg*m^2
J = I_S_Sstar(1);  % kg*m^2
p = H/I  % precession rate, rad/s
s = (I-J)/I * N_omega_S(1)  % spin rate, rad/s
phi = acosd(J*s/((I-J)*p))  % nutation angle, deg

% calculate angles at t = 0.2 s
beta_t02 = rad2deg(p*0.2)  % precession angle, deg
psi_t02 = rad2deg(s*0.2)  % spin angle, deg

% calculate angles at t = 1 s
beta_t10 = rad2deg(p)  % precession angle, deg
psi_t10 = rad2deg(s)  % spin angle, deg

fprintf('\nPart c:\n\n');

% calculate Euler parameters at t = 0.2 s
psi = psi_t02; beta = beta_t02;
N_curlE_S = zeros(1,4);  % vector component in c_hat frame
N_curlE_S(1) = cosd(phi)*sind(beta/2)*cosd(psi/2) + cosd(beta/2)*sind(psi/2)
N_curlE_S(2) = sind(phi)*sind(psi/2)*sind(beta/2)
N_curlE_S(3) = -sind(phi)*sind(beta/2)*cosd(psi/2)
N_curlE_S(4) = cosd(beta/2)*cosd(psi/2) - cosd(phi)*sind(beta/2)*sind(psi/2)
N_curlE_S_check = sum(N_curlE_S.^2)  % should be 1 (it is!)
C_C_S = [1 0 0; 0, cosd(psi), -sind(psi); 0, sind(psi), cosd(psi)]
N_curlE_S(1:3) = N_curlE_S(1:3)*C_C_S

% calculate Euler parameters at t = 1.0 s
psi = psi_t10; beta = beta_t10;
N_curlE_S = zeros(1,4);  % vector component in c_hat frame
N_curlE_S(1) = cosd(phi)*sind(beta/2)*cosd(psi/2) + cosd(beta/2)*sind(psi/2)
N_curlE_S(2) = sind(phi)*sind(psi/2)*sind(beta/2)
N_curlE_S(3) = -sind(phi)*sind(beta/2)*cosd(psi/2)
N_curlE_S(4) = cosd(beta/2)*cosd(psi/2) - cosd(phi)*sind(beta/2)*sind(psi/2)
N_curlE_S_check = sum(N_curlE_S.^2)  % should be 1 (it is!)
C_C_S = [1 0 0; 0, cosd(psi), -sind(psi); 0, sind(psi), cosd(psi)]
N_curlE_S(1:3) = N_curlE_S(1:3)*C_C_S

% Part e: plot inertia ellipsoid projection again,
%         for drawing space/body cones
figure();
plot(x1, x3, '-b'); hold on; plot(x1, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 08, Problem 1e: Inertia Ellipsoid Projection (s1-s3)')
xlabel('s_1'); ylabel('s_3'); grid on; hold on;
axis square; axis([-1 1 -1 1] * scale);

%% Problem 3
fprintf('\nProblem 3:\n');

fprintf('\nPart a:\n\n');

I = [500 252 252]  % principal MOI's, kg*m^2

% calculate inertia ellipsoid semi-diameters on corresponding axes
alpha = I.^(-1/2)

% plot inertia ellipsoid projections
scale = [-1 1 -1 1] * 0.065;

% s1-s3 projection
x1 = linspace(-alpha(1),alpha(1),1000);
x3 = sqrt(1/I(3) * (1 - I(1)*x1.^2));
figure();
plot(x1, x3, '-b'); hold on; plot(x1, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 08, Problem 2a: Inertia Ellipsoid Projection (s1-s3)')
xlabel('s_1'); ylabel('s_3'); grid on;
axis square; axis(scale)

% s1-s2 projection
x2 = sqrt(1/I(2) * (1 - I(1)*x1.^2));
figure();
plot(x1, x2, '-b'); hold on; plot(x1, -x2, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 08, Problem 2a: Inertia Ellipsoid Projection (s1-s2)')
xlabel('s_1'); ylabel('s_2'); grid on;
axis square; axis(scale)

% s2-s3 projection
x2 = linspace(-alpha(2),alpha(2),1000);
x3 = sqrt(1/I(3) * (1 - I(2)*x2.^2));
figure();
plot(x2, x3, '-b'); hold on; plot(x2, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 08, Problem 2a: Inertia Ellipsoid Projection (s2-s3)')
xlabel('s_2'); ylabel('s_3'); grid on;
axis square; axis(scale)

fprintf('\nPart b:\n\n');

N_H_S_Sstar = [I(1)*N_omega_S(1), 0, I(3)*N_omega_S(3)]
  % angular momentum vector in s_hat frame, kg*m^2/s
phi = atand(N_H_S_Sstar(3)/N_H_S_Sstar(1))  % nutation angle, deg

fprintf('\nPart c:\n\n');

p = norm(N_H_S_Sstar)/I(2)  % precession rate, rad/s
s = (I(2) - I(1))/I(2) * N_omega_S(1)  % spin rate, rad/s

%% Problem 3
fprintf('\nProblem 3:\n\n');

% numerically simulate dynamic and kinematic differential equations
% simultaneously
% NOTE: this code gets results for parts a, b, and d together
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
  % error tolerance is set to 10^(-7)
[t, y] = ode45('hw08_numint_all', t_span, y_0, options);

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

n = size(y,1);  % number of output stats
beta = zeros(n,1);  % precession angle, deg
sigma = zeros(n,1);  % nutation angle, deg
beta(1) = 0; sigma(1) = 0;  % initial conditions
                            % (avoid divide by zero in for loop)

k = 1;
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
end

% create plots

% plot angular velocities
figure();
plot(t, w1, '-b'); hold on; plot(t, w2, '-r');
plot(t, w3, '-g');
title('Jordan Mayer - AAE 590/440, HW 08, Problem 3a: Zero-Torque Results');
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)'); grid on;
ylim([-2.5 2.5]);
legend('s-hat_1 component', 's-hat_2 component', 's-hat_3 component');

% plot precession and nutation angles (beta, sigma)
figure();
plot(t, beta, '-b'); hold on; plot(t, sigma, '-r');
title('Jordan Mayer - AAE 590/440, HW 08, Problem 3b: Zero-Torque Results');
xlabel('Time (s)'); ylabel('Angle (deg)'); grid on;
legend('precession (beta)', 'nutation (sigma)');

% plot direction cosine components
figure();
plot(C21, C31, '-b');
title('Jordan Mayer - AAE 590/440, HW 08, Problem 3d: Zero-Torque Results');
xlabel('C_{21}'); ylabel('C_{31}'); grid on;
axis square; axis([-1 1 -1 1]);

% get results for times 0.2 s and 1.0 s
ind_t02 = find(abs(t - 0.2) < 0.005);  % index for t = 0.2 s
C21_t02 = C21(ind_t02)
C31_t02 = C31(ind_t02)
beta_t02 = beta(ind_t02)
ind_t10 = find(abs(t - 1.0) < 0.005);  % index for t = 1.0 s
C21_t10 = C21(ind_t10)
C31_t10 = C31(ind_t10)
beta_t10 = beta(ind_t10)