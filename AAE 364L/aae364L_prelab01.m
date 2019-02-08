%%%%%
% Jordan Mayer
% AAE 364L
% Prelab 01
% January 16, 2019
%%%%%

clear all; close all; format compact;

% Part (ii)

m = 1.0731;  % kg
c = 13.124;  % kg/s
k = 10;  % V/m
gamma = 1.72;  % s.A/m

H = tf([gamma], [m c gamma*k 0])  % open-loop transfer function
rlocus(H)
title('Jordan Mayer, Prelab 01 - Root Locus of H(s)');

options = simset('SrcWorkspace', 'current');
sim('aae364L_prelab01_sys', [], options);

S = stepinfo(y,t)
e = y(size(y,1)) - 0.5