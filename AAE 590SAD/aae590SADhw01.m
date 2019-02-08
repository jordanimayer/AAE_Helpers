%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
% HW 01
% January 18, 2019
%
% Perform basic vector/matrix arithmetic to obtain and verify answers.
%%%%%

clear all; close all; format compact;

%%%%% Problem 3 %%%%%

fprintf('\n----- Problem 3 -----\n');

%%% Part i %%%

fprintf('\n--- Part i ---\n\n');

% express all vectors as column vectors in y_hat frame
% ignore d since this will cancel out when converting to unit vectors
b1 = [-6 4 0]'
b_hat1 = b1/norm(b1)
v = [0 -4 -3]'
b2 = cross(v, b1)
b_hat2 = b2/norm(b2)
b_hat3 = cross(b_hat1, b_hat2)

L = [-0.8321 0.3714 -0.4120; 0.5547 0.5571 -0.6180; 0 -0.7428 -0.6695]
C = L.'
orth_test = L*C

%%% Part ii %%%

fprintf('\n--- Part ii ---\n\n');

% still express as column vector in y_hat frame
q = [2 4/3 -1]'
q_hat_y = q/norm(q)

% now, in b_hat frame...
q_hat_b = C * q_hat_y

% no sign convention yet, for now just saw all these angles are positive
ang1 = acosd(-0.3551)
ang2 = acosd(0.8560)
ang3 = acosd(-0.3759)