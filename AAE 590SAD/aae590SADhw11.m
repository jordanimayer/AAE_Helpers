%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 11
%
% Math! Numerical simulation! Yeah!
%%%%%

%% Preliminary setup

close all; clear all; format compact; format short; rehash toolbox;
cd('W:\Courses Spring 2019\AAE 590SAD\AAE 590SAD Homework\HW 11');
addpath('../helpers');

%% Problem 1b
fprintf('\nProblem 1b:\n');

I_A = 361.25;  % kg-met^2
I_B = 861.25;  % kg-met^2
I_C = 1062.5;  % kg-met^2
I = zeros(1,3);  % will hold moments of inertia w.r.t. a-frame

fprintf('\nPart i:\n\n');

I = [I_A, I_C, I_B];
Ki = get_K_from_I(I)

fprintf('\nPart ii:\n\n');

I = [I_C, I_A, I_B];
Kii = get_K_from_I(I)

fprintf('\nPart iii:\n\n');

I = [I_B, I_A, I_C];
Kiii = get_K_from_I(I)

fprintf('\nPart iv:\n\n');

I = [I_B, I_C, I_A];
Kiv = get_K_from_I(I)

%% Problem 3
clear all; close all;

fprintf('\nProblem 3:\n');

% set up constants
w1_0 = 1.0;
  % initial nondimensionalized angular velocity parallel to orbit normal
  % (B in N)
q_0 = [0.0, 0.0, 0.0, 1.0];
  % initial Euler parameters (A to B)
v_span = [0, 10];  % span of revolutions to simulate for
region_list = [1, 4, 6, 3];  % orientation regions to simulate
w23_0_list = [0, 0.08, 0.16];
  % initial nondimensionalized angular velocity parallel to orbit radius
  % (w2) and parallel to orbit tangent (w3)
  % (B in N) (all three conditions)

% simulate, plot angles, and get eigenvalues
for region = region_list
    fprintf(['\nRegion ', num2str(region), ':\n\n']);
    global K;  % defined in hw11_sim_and_plot for each region
    
    for w23_0 = w23_0_list
        w_0 = [w1_0, w23_0, w23_0];
        hw11_sim_and_plot(v_span, region, w_0, q_0);
    end
    
    fprintf(['K = ', num2str(K), '\n\n']);
    
    mu1 = 0;
    mu2 = sqrt(3*K(1));
    mu3 = -mu2;
    b = (1 - K(2)*K(3) + 3*K(3))/2;
    c = -4*K(2)*K(3);
    mu4 = (-b + sqrt(b^2 - c))^(1/2);
    mu5 = -mu4;
    mu6 = (-b - sqrt(b^2 - c))^(1/2);
    mu7 = -mu6;
    mu = [mu1, mu2, mu3, mu4, mu5, mu6, mu7].'
end

%% Problem 4, quick part
clear all; close all;

fprintf('\nProblem 4:\n\n');

% find moments of inertia for total gyrostat for each region
global I_reg1 I_reg4 I_reg6;

I_A = 361.25;  % kg-met^2
I_B = 861.25;  % kg-met^2
I_C = 1062.5;  % kg-met^2
I_rot = [82.50, 96.25, 96.25];  % rotor moments of inertia, kg-met^2

I_reg1 = [I_C, I_A, I_B] + I_rot
I_reg4 = [I_A, I_C, I_B] + I_rot
I_reg6 = [I_B, I_A, I_C] + I_rot

%% Problem 4, simulations/plots
close all;

% set up constants
w1_0 = 1.0;
  % initial nondimensionalized angular velocity parallel to orbit normal
  % (B in N)
q_0 = [0.0, 0.0, 0.0, 1.0];
  % initial Euler parameters (A to B)
v_span = [0, 10];  % span of revolutions to simulate for
region_list = [1, 4, 6];  % orientation regions to simulate
w23_0_list = [0, 0.08, 0.16];
  % initial nondimensionalized angular velocity parallel to orbit radius
  % (w2) and parallel to orbit tangent (w3)
  % (B in N) (all three conditions)
global n_rotor;
n_rotor = 50;   % rotor spin rate, nondimensionalized with respect to
                  % orbit rate
 
% simulate and plot angles
for region = region_list
    for w23_0 = w23_0_list
        w_0 = [w1_0, w23_0, w23_0];
        hw11p4_sim_and_plot(v_span, region, w_0, q_0);
    end
end