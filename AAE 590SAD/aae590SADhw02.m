%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 02
%
% Perform basic matrix arithmetic and call functions to evaluate dyads.
%%%%%

clear all; close all; rehash toolbox; format compact;

%% Problem 1

fprintf('\nProblem 1:\n\n');

nCb = [.8925 .1574 -.4226; -.2754 .9323 -.2346; .3571 .3258 .8754];
  % direction cosine matrix from n-frame to b-frame
orth_check = nCb * nCb.'  % should be identity matrix (it is!)

%% Problem 2

fprintf('\nProblem 2:\n\n');

Ib = [100 0 -100; 0 500 0; -100 0 400];
  % inertia dyadic in b-frame, matrix form (kg/m^2)
In_mat = nCb * Ib * nCb.'  % inertia dyadic in n-frame, matrix form
                           % as evaluated using similarity transformation
                           % (kg*m^2)
In_dyad = convert_dyadic_dircos(nCb, Ib)
  % inertia dyadic in n-frame, matrix form
  % as evaluated using direct dyadic
  % (kg*m^2)
  
w = [2 0 -1];  % angular velocity vector of n-frame relative to b-frame,
               % expressed in n-frame (rad/s)
H = dyadic_vec_dot(In_dyad, w)
  % angular momentum vector expressed in n-frame, matrix form (kg*m^2/s)

T = dot(w/2, H)  % rotational kinetic energy (kg*m^2/s^2)

%% Problem 3

fprintf('\nProblem 3:\n\n');

L = [2 -2 1]  % rotation axis, in a-frame (or b-frame)
lambda = L/norm(L)  % rotation axis unit vector, in a-frame (or b-frame)
theta = 60;  % rotation angle, deg

aCb = dircos_lambdatheta(lambda, theta)
orth_check = aCb * aCb.'  % should be identity matrix (it is!)
aRb = srd_lambdatheta(lambda, theta)

k_a = [3 1 0];  % original vector, in a-frame
k_b_a = rotvec_srt(k_a, lambda, theta)  % rotated vector, in a-frame