%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 06
%
% Perform basic arithmetic to obtain and verify results.
%%%%%

%% Preliminary setup
close all; clear all; format compact; format long; rehash toolbox;
addpath('../helpers');
G = 6.67e-11;  % m^3/(kg*s)

%% Problem 1
fprintf('\nProblem 1:\n\n');

% all masses in terms of m (in kg)
m1 = 3; m2 = 4; m3 = 2; m4 = 5; mP = 5;

% all distances in terms of k (in meters), all vectors in the n_hat frame
r1 = [-2 -4 0]; r2 = [1 -2 0]; r3 = [2 -6 0]; r4 = [3 5 0]; rP = [-2 2 0];

fprintf('\nPart a:\n\n');

% get c.m. location
rCM = (m1*r1 + m2*r2 + m3*r3 + m4*r4)/(m1+m2+m3+m4)

fprintf('\nPart b:\n\n');

% get relative positions (in terms of k) and gravity forces (in terms of
% (m/k)^2
mS = [m1 m2 m3 m4]; rS = [r1; r2; r3; r4];
for i=1:4
    p(i,:) = rS(i,:) - rP
    norm_p(i) = norm(p(i,:))
    F(i,:) = -G*mP*mS(i)*p(i,:)*((norm_p(i))^2)^(-3/2) * 10^12
      % pico-Newtons
end

% get total gravity force
Fg = F(1,:) + F(2,:) + F(3,:) + F(4,:)  % pico-Newtons
norm_Fg = norm(Fg)
Fg_hat = Fg/norm_Fg

fprintf('\nPart c:\n\n');

% get distance from P' to c.g. (in terms of k)
RP = (G * mP * sum(mS)/(norm_Fg*10^(-12)))^(1/2)
rPcg = -RP*Fg_hat  % P' to c.g. in terms of k, vector

fprintf('\nPart d.i:\n\n');

% get vector from c.m. to P', in terms of k
rCM_P = rP - rCM

% get moment about c.m. due to Fg acting at P' in pico-Newton-meters, in
% terms of m^2/k
M_CM = cross(rCM_P, Fg)

fprintf('\nPart d.ii:\n\n');

% get vector from c.m. to c.g., in terms of k
rCM_CG = rCM_P + rPcg

% get moment about c.m. due to Fg acting at c.g. in pico-Newton-meters, in
% terms of m^2/k
M_CM = cross(rCM_CG, Fg)

%% Problem 2
fprintf('\nProblem 2:\n\n');

% all distances in terms of L, measured from P' in n_hat frame as defined
% on handwritten sketch (all points defined in handwritten work)
r1 = [1 0 0];
r2 = [1+cosd(60), sind(60) 0];
r3 = [1+cosd(60), -sind(60) 0];

% all masses in terms of m
mP = 2; m1 = 1/2; m2 = 1; m3 = 1;

fprintf('\nPart a:\n\n');

% get c.m. location (in terms of L)
rCM = (m1*r1 + m2*r2 + m3*r3)/(m1+m2+m3)

% get relative positions (in terms of L) and gravity forces (in terms of
% (m/L)^2
mB = [m1 m2 m3]; rB = [r1; r2; r3];
for i=1:3
    p(i,:) = rB(i,:)
    norm_p(i) = norm(p(i,:))
    F(i,:) = -G*mP*mB(i)*p(i,:)*((norm_p(i))^2)^(-3/2) * 10^12
      % pico-Newtons
end

% get total gravity force (in terms of (m/L)^2), pico-Newtons,
% vector in n_hat frame
Fg = F(1,:) + F(2,:) + F(3,:)
norm_Fg = norm(Fg)
Fg_hat = Fg/norm_Fg

% get distance from P' to c.g. (in terms of L)
RP = (G * mP * sum(mB)/(norm_Fg*10^(-12)))^(1/2)
rPcg = -RP*Fg_hat  % P' to c.g. in terms of L, vector in n_hat frame

fprintf('\nPart d:\n\n');

% get new c.m. location (in terms of L
rCM = (m1*r1 + m2*r2)/(m1+m2)

% get new total gravity force (in terms of (m/L)^2), pic-Newtons,
% vector in n_hat frame
Fg = F(1,:) + F(2,:)
norm_Fg = norm(Fg)
Fg_hat = Fg/norm_Fg

% get distance from P' to c.g. (in terms of L)
RP = (G*mP*(m1+m2)/(norm_Fg*10^(-12)))^(1/2)
rPcg = -RP*Fg_hat  % P' to c.g. in terms of L, vector in n_hat frame

% get vector from c.m. to c.g., in terms of k
rCM_CG = rPcg - rCM

% get moment about c.m. due to Fg acting at c.g. in pico-Newton-meters, in
% terms of m^2/k
M_CM = cross(rCM_CG, Fg)

%% Problem 3
fprintf('\nProblem 3:\n\n');

fprintf('\nPart b:\n\n');

% moments of inertia for each piece (about B*), in terms of m*L^2
I_1 = [1/6 1/32 1/6];  % [I_11, I_22, I_33]
I_2 = [1/4 7/192 .2866];
I_3 = I_2;

% full body inertia matrix (about B*), in terms of m*L^2
I = diag(I_1 + I_2 + I_3)

fprintf('\nPart c:\n\n');

mP = 5.97e24;  % kg
R = 6488*1000; L = 70*1000;  % meters

tr_I = trace(I)  % trace of I, in terms of m*L^2
f_2_curly = 3/2*(1.511 - 5*.667) + 3*.667
  % term inside curly braces of f_2 equation, in terms of m*L^2,
  % in b_hat1 direction
f_2 = f_2_curly*L^2/R^2
F = -G*mP/R^2 * (1 + f_2)
  % total gravity force, in terms of m, in a_hat1 direction, Newtons

% Part d
Rcg = (G*mP/abs(F))^(1/2)  % distance from Earth to c.g., meters
Rcg = Rcg/1000  % km

% Part e
f_2_curly_1 = 0.5*(.6667*(1-3*.7071) + .1042*(1-3*(-.7071)))
  % a_hat1 term of f_2 term in curly braces
f_2_curly_2 = .6667*.7071^2 + .1042*.7071^2*(-1)
  % a_hat2 term of f_2 term in curly braces
f_2_1 = 3*(70^2/6488^2)*f_2_curly_1  % a_hat1 term of f_2
f_2_2 = 3*(70^2/6488^2)*f_2_curly_2  % a_hat2 term of f_2
F = -G*mP/R^2 * [1+f_2_1, f_2_2, 0]
  % total gravity force, in terms of m, in a_hat vector frame, Newtons
norm_F = norm(F)  % in terms of m, Newtons
Rcg = (G*mP/norm_F)^(1/2)  % distance from Earth to c.g., meters
Rcg = Rcg/1000  % km
F_hat = F/norm_F
  % unit vector in gravity force direction, in a_hat vector frame
rP_cg = -Rcg*F_hat
  % vector from Earth to c.g., in a_hat vector frame, km
rP_cm = [6387.988, 12.37, 0]
  % vector from Earth to c.m., in a_hat vector frame, km
Mcm = cross(rP_cm, F)
  % moment about c.m. due to gravity force, in terms of m, Newton-km