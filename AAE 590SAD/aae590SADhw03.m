%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 03
% February 1, 2019
%
% Perform matrix arithmetic to obtain and verify results.
%%%%%

close all; clear all; format compact; rehash toolbox;
addpath('../helpers');

%% Problem 1
fprintf('\nProblem 1:\n\n');

N_C_B0 = [0.9363 -0.2896 0.1987; 0.3130 0.9447 -0.0981; -0.1593 0.1540 0.9751];
B0_lambda_B_n = [2 3 6]/7  % in n-frame
B0_theta_B = 30;  % deg
B0_lambda_B_b0 = B0_lambda_B_n * N_C_B0  % in b0-frame
B0_C_B = dircos_lambdatheta(B0_lambda_B_b0, B0_theta_B)
N_C_B = N_C_B0 * B0_C_B

N_q_B = quat_from_dircos(N_C_B)  % vector part in either n- or b-frame

[N_lambda_B, N_theta_B] = ser_from_quat(N_q_B)
  % lambda in either n- or b-frame
  % theta in deg
check = N_C_B - dircos_lambdatheta(N_lambda_B, N_theta_B)
  % should be near zero (it is!)
  
%% Problem 2
fprintf('\nProblem 2:\n\n');

q = [cosd(60) sind(60) 0]  % o-frame
O_lambda_Sprime = q;
O_theta_Sprime = 135;  % deg
O_C_Sprime = dircos_lambdatheta(O_lambda_Sprime, O_theta_Sprime)

Sprime_lambda_S = [0 1 0];  % s-frame
Sprime_theta_S = 30;  % deg
Sprime_C_S = dircos_lambdatheta(Sprime_lambda_S, Sprime_theta_S)

O_C_S = O_C_Sprime * Sprime_C_S

O_q_Sprime = quat_from_dircos(O_C_Sprime)
  % vector part in O-frame or Sprime-frame
Sprime_q_S = quat_from_dircos(Sprime_C_S)
  % vector part in Sprime-frame or S-frame
O_q_S = quat_from_dircos(O_C_S)
  % vector part in O-frame or S-frame
  
O_q_S = quat_from_qsequence(O_q_Sprime, Sprime_q_S)
  % vector part in Sprime-frame
O_q_S(1:3) = O_q_S(1:3) * Sprime_C_S
  % vector part now in S-frame or O-frame
  
[O_lambda_S, O_theta_S] = ser_from_quat(O_q_S)
  % lambda in either O- or S-frame
  % theta in deg
check = O_C_S - dircos_lambdatheta(O_lambda_S, O_theta_S)
  % should be near zero (it is!)

%% Problem 3
fprintf('\nProblem 3:\n\n');

phi212 = [200 -75 30]  % rotation angles, 2-1-2, deg
N_C_B = dircos_b212(phi212)

N_q_B = quat_from_dircos(N_C_B)
  % vector component in N- or B-frame

[N_lambda_B, N_theta_B] = ser_from_quat(N_q_B)
  % lambda in N- or B-frame
  % theta in deg
  
theta2_b321 = asind(-N_C_B(3,1))  % deg, quadrant by convention
theta1_b321_1 = asind(N_C_B(2,1)/cosd(theta2_b321))  % deg, option 1
theta1_b321_2 = 180 - theta1_b321_1 - 360  % deg, option 2, -360 to match
theta1_b321_3 = acosd(N_C_B(1,1)/cosd(theta2_b321))  % deg, option 3
theta1_b321_4 = -theta1_b321_3  % deg, option 4
theta1_b321 = theta1_b321_2  % deg, chosen for consistency
theta3_b321_1 = asind(N_C_B(3,2)/cosd(theta2_b321))  % deg, option 1
theta3_b321_2 = 180 - theta3_b321_1  % deg, option 2
theta3_b321_3 = acosd(N_C_B(3,3)/cosd(theta2_b321))  % deg, option 3
theta3_b321_4 = -theta3_b321_3  % deg, option 4
theta3_b321 = theta3_b321_2  % deg, chosen for consistency
phi1_b321 = theta1_b321  % deg
phi2_b321 = theta2_b321  % deg
phi3_b321 = theta3_b321  % deg

phi2_s212 = acosd(N_C_B(2,2))  % deg, quadrant by convention
phi1_s212_1 = asind(N_C_B(2,1)/sind(phi2_s212))  % deg, option 1
phi1_s212_2 = 180 - phi1_s212_1 - 360  % deg, option 2, -360 to match
phi1_s212_3 = acosd(-N_C_B(2,3)/sind(phi2_s212))  % deg, option 3
phi1_s212_4 = -phi1_s212_3  % deg, option 4
phi1_s212 = phi1_s212_2  % deg, chosen for consistency
phi3_s212_1 = asind(N_C_B(1,2)/sind(phi2_s212))  % deg, option 1
phi3_s212_2 = 180 - phi3_s212_1  % deg, option 2
phi3_s212_3 = acosd(N_C_B(3,2)/sind(phi2_s212))  % deg, option 3
phi3_s212_4 = -phi3_s212_3  % deg, option 4
phi3_s212 = phi3_s212_1  % deg, chosen for consistency