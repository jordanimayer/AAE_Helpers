%%%%%
% Jordan Mayer
% AAE 364L
% Prelab 05
%%%%%

%% Preliminary setup

close all; clear all; format compact; format long; rehash toolbox;
addpath('W:\Courses Spring 2019\AAE 364L\AAE 364L Labs\Lab 05\Lab 05 Files');

%% Part (i)

run('setup_heli_PID_parameters');

l = 0.185;
k_p = l*K_pp/0.185/J_eq_p
a_p = B_p/J_eq_p
b_p = m_heli*g*h/J_eq_p

lambda_p = [-1 -3+3i, -3-3i]
p_p = poly(lambda_p)
alpha_p = p_p(1)
beta_p = p_p(2)
gamma_p = p_p(3)

gamma_dp = (alpha_p - a_p)/k_p
gamma_pp = (beta_p - b_p)/k_p
gamma_ip = gamma_p/k_p

% Part (ii)

k_y = l*K_yy/0.185/J_eq_y
a_y = B_y/J_eq_y
b_y = 0

w_n = 1;  % rad/s
zeta = 0.7;
lambda = -zeta*w_n + (w_n*sqrt(1-zeta^2))*1i
z = k_y/(lambda^2 + a_y*lambda + b_y)

gamma_iy = 1;

M = [lambda^2 * z, lambda*z; conj(lambda)^2 * conj(z), conj(lambda) * conj(z)];
N = [lambda + gamma_iy*z; conj(lambda) + gamma_iy*conj(z)];
gammas_y = -real(inv(M)*N);
gamma_dy = gammas_y(1)
gamma_py = gammas_y(2)

% Part (iii)
gamma_pp = 39;
gamma_py = 39;
gamma_dp = 39;
gamma_dy = 39;
gamma_ip = 4;
gamma_iy = 4;
Ki = [gamma_pp, 0, gamma_dp, 0, gamma_ip, 0; 0, gamma_py, 0, gamma_dy, 0, gamma_iy]