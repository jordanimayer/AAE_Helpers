%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 07
%
% Perform basic arithmetic to obtain and verify results.
% Plot inertia ellipsoids both as 2-D projections and in 3-D.
%%%%%

%% Preliminary setup
close all; clear all; format compact; rehash toolbox;
addpath('../helpers');

%% Problem 1

% Part a

I = [150 500 500];  % principal moments of inertia, kg*m^2

% calculate semi-diameters on corresponding axes
alpha = 1./sqrt(I)
scale = [-1 1 -1 1] * 0.09;

% plot projections

% s1-s3
x1 = linspace(-alpha(1),alpha(1),1000);
x3 = sqrt(1/I(3) * (1 - I(1)*x1.^2));
figure();
plot(x1, x3, '-b'); hold on; plot(x1, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid Projection (s1-s3)')
xlabel('s_1'); ylabel('s_3'); grid on;
axis square; axis(scale)

% s1-s2
x2 = sqrt(1/I(2) * (1 - I(1)*x1.^2));
figure();
plot(x1, x2, '-b'); hold on; plot(x1, -x2, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid Projection (s1-s2)')
xlabel('s_1'); ylabel('s_2'); grid on;
axis square; axis(scale)

% s2-s3
x2 = linspace(-alpha(2),alpha(2),1000);
x3 = sqrt(1/I(3) * (1 - I(2)*x2.^2));
figure();
plot(x2, x3, '-b'); hold on; plot(x2, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid Projection (s2-s3)')
xlabel('s_2'); ylabel('s_3'); grid on;
axis square; axis(scale)

% 3-D plot
[x,y,z] = ellipsoid(0, 0, 0, alpha(1), alpha(2), alpha(3), 40);
figure();
surf(x,y,z);
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid (3-D)')
xlabel('s_1'); ylabel('s_2'); zlabel('s_3');
axis square; axis([scale scale(1:2)])

% Part b

I = [600 375 375];  % principal moments of inertia, kg*m^2

% calculate semi-diameters on corresponding axes
alpha = 1./sqrt(I)
scale = [-1 1 -1 1] * 0.06;

% plot projections

% s1-s3
x1 = linspace(-alpha(1),alpha(1),1000);
x3 = sqrt(1/I(3) * (1 - I(1)*x1.^2));
figure();
plot(x1, x3, '-b'); hold on; plot(x1, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid Projection (s1-s3)')
xlabel('s_1'); ylabel('s_3'); grid on;
axis square; axis(scale)

% s1-s2
x2 = sqrt(1/I(2) * (1 - I(1)*x1.^2));
figure();
plot(x1, x2, '-b'); hold on; plot(x1, -x2, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid Projection (s1-s2)')
xlabel('s_1'); ylabel('s_2'); grid on;
axis square; axis(scale)

% s2-s3
x2 = linspace(-alpha(2),alpha(2),1000);
x3 = sqrt(1/I(3) * (1 - I(2)*x2.^2));
figure();
plot(x2, x3, '-b'); hold on; plot(x2, -x3, '-b');  % include negative sqrts
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid Projection (s2-s3)')
xlabel('s_2'); ylabel('s_3'); grid on;
axis square; axis(scale)

% 3-D plot
[x,y,z] = ellipsoid(0, 0, 0, alpha(1), alpha(2), alpha(3), 40);
figure();
surf(x,y,z);
title('Jordan Mayer - AAE 590/440, HW 07, Problem 1a: Inertia Ellipsoid (3-D)')
xlabel('s_1'); ylabel('s_2'); zlabel('s_3');
axis square; axis([scale scale(1:2)])

%% Problem 2

% Part c

G = 6.67e-11  % gravitational constant, m^3/(kg*s^2)

% relative position vectors of particles from P, in terms of L, expressed
% in a_hat frame
p1 = [1.5 .5 -.5]
p2 = [1.5 .5 .5]
p3 = [1.5 -.5 -.5]
p4 = [1.5 -.5 .5]
p5 = [2.5 .5 -.5]
p6 = [2.5 .5 .5]
p7 = [2.5 -.5 -.5]
p8 = [2.5 -.5 .5]
p_all = [p1; p2; p3; p4; p5; p6; p7; p8]

% gravitational force on each particle due to P, in terms of m'*m/L^2,
% Newtons, expressed in a_hat frame
F_all = zeros(8,3)

% total gravitational force, in terms of m'*m/L^2, Newtons, expressed in
% a_hat frame
F = zeros(1,3)

% calculate gravitational forces from each particle and sum to get total
% force
for k=1:8
    p_k = p_all(k,:)
    p_k_norm = norm(p_k)
    F_all(k,:) = -(G/p_k_norm^3)*p_k
    F = F + F_all(k,:)
end

% Part d

% magnitude of total gravitational force, in terms of m'*m/L^2, Newtons
F_norm = norm(F)

% unit vector in direction of total gravitational force, expressed in a_hat
% frame
F_hat = F/F_norm

% distance from P to c.g., in terms of L
R_cg = sqrt(G*8/F_norm)

% vector from P to c.g., in terms of L, in a_hat frame
r_P_cg = -R_cg*F_hat

% vector from c.m. to c.g., in terms of L, in a_hat frame
r_cm_cg = r_P_cg - [2 0 0]

% distance from c.m. to c.g., in terms of L
r_cm_cg_norm = norm(r_cm_cg)