%%%%%
% Jordan Mayer
% AAE 364L
% Prelab 06
%%%%%

%% Preliminary setup

close all; clear all; format compact; format short; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 06\Lab 06 Files');

%% Finding gains using LQR
close all;

run('setup_heli_LQR_parameters');

% form Ai, Bi (they're already in the setup file!)
Ai = a; Bi = b;

% use LQR to find Ki for both helicopter models


% simulate helicopter using error feedback and adjust weights/gains as
% needed
Q = diag([9 9 1 1 9 9]); R = diag([0.01, 0.02]);
Ki = lqr(Ai, Bi, Q, R)
sim('heli_model');

% simulate helicopter using state feedback and adjust weights/gains as
% needed
Q = diag([7 7 3 3 9 9]); R = diag([0.01, 0.05]);
Ki = lqr(Ai, Bi, Q, R)
sim('heli_model2');