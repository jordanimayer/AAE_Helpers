%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 04
%
% Perform basic arithmetic to obtain and verify results.
%%%%%

close all; clear all; format compact; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 590SAD\AAE 590SAD Homework\helpers');

%% Problem 2

w = [1 -3 -1.5];  % rad/s, in n-frame
q = [.25 .5 -.75 0];  % Euler parameter, vector component in n-frame
q(4) = sqrt(1 - q(1)^2 - q(2)^2 - q(3)^2)
temp1 = q(4)*w
temp2 = cross(q(1:3), w)
q_dot = zeros(1,4); % time-derivative of Euler parameter relative to
                    % B (vector component in n-frame

%% Part (i)

C = dircos_quat(q)  % direction cosine matrix from N to B

% all angles in deg
beta_2 = acosd(C(1,1))
s2 = sind(beta_2)
beta_1 = [asind(C(1,2)/s2), 180 - asind(C(1,2)/s2), acosd(C(1,3)/s2), -acosd(C(1,3)/s2)]
beta_3 = [asind(C(2,1)/s2), 180 - asind(C(2,1)/s2), acosd(C(3,1)/s2), -acosd(C(3,1)/s2)]
beta_1 = beta_1(2)
beta_3 = beta_3(4)

% all rates in rad/s
beta_dot_3 = (w(2)*sind(beta_1) + w(3)*cosd(beta_1))/sind(beta_2)
beta_dot_2 = (w(2) - beta_dot_3*sind(beta_1)*sind(beta_2))/cosd(beta_1)
beta_dot_1 = w(1) - beta_dot_3*cosd(beta_2)

%% Part (ii)

C = C  % same direction cosine matrix

% all angles in deg
beta_2 = asind(C(1,3))
c2 = cosd(beta_2)
beta_3 = [acosd(C(1,1)/c2), -acosd(C(1,1)/c2), asind(-C(1,2)/c2), 180-asind(-C(1,2)/c2)]
beta_1 = [asind(-C(2,3)/c2), 180-asind(-C(2,3)/c2), acosd(C(3,3)/c2), -acosd(C(3,3)/c2)]
beta_3 = beta_3(1)
beta_1 = beta_1(1)

% all rates in rad/s
beta_dot_2 = w(2)*cosd(beta_3) + w(1)*sind(beta_3)
beta_dot_1 = (w(1) - beta_dot_2*sind(beta_3))/(cosd(beta_2)*cosd(beta_3))
beta_dot_3 = w(3) - beta_dot_1*sind(beta_2)

%% Problem 3

% all angles in deg
nu1 = -250; nu2 = 60; nu3 = 90;
% all rates in rad/s
nu_dot1 = 0; nu_dot2 = 1; nu_dot3 = -2;

%% Part b

s1 = sind(nu1)
s2 = sind(nu2)
s3 = sind(nu3)
c1 = cosd(nu1)
c2 = cosd(nu2)
c3 = cosd(nu3)

C = zeros(3,3);
C(1,1) = s1*s2*s3 + c3*c1;
C(1,2) = c1*s2*s3 - c3*s1;
C(1,3) = c2*s3;
C(2,1) = s1*c2;
C(2,2) = c1*c2;
C(2,3) = -s2;
C(3,1) = s1*s2*c3 - s3*c1;
C(3,2) = c1*s2*c3 + s3*s1;
C(3,3) = c2*c3

beta = 90
gamma = [acosd(C(2,3)), -acosd(C(2,3)), asind(C(1,3)), 180-asind(C(1,3))]
phi = [acosd(C(3,2)), -acosd(C(3,2)), asind(C(3,1)), 180-asind(C(3,1))]
gamma = gamma(1)
phi = phi(1)

% all rates in rad/s
w = zeros(1,3);
w(1) = nu_dot2*c1 + nu_dot3*s1*c2;
w(2) = -nu_dot2*s1 + nu_dot3*c1*c2;
w(3) = nu_dot1 - nu_dot3*s2

% now use the body 1-3-1 angles
s2 = sind(beta)
s3 = sind(phi)
c2 = cosd(beta)
c3 = cosd(phi)

% all rates in rad/s
gamma_dot = (-w(2)*c3 + w(3)*s3)/s2
beta_dot = w(2)*s3 + w(3)*c3
phi_dot = w(1) + (w(2)*c3 - w(3)*s3)*c2/s2

q = quat_from_dircos(C)