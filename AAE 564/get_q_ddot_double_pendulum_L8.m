%%%%%
% get_q_ddot_double_pendulum
%
% Get q_ddot for double-pendulum problem, using linearized equation for
% L8 conditions (see HW 06)
%
% Inputs:
%   u
%   q: [y, theta1, theta2]'
%   q_dot: [y_dot, theta1_dot, theta2_dot]'
%
% Outputs:
%   q_ddot: [y_ddot, theta1_ddot, theta2_ddot]
%%%%%

function q_ddot = get_q_ddot_double_pendulum_L8(u, q, q_dot)

x = [q; q_dot];  % state vector

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 -0.5 -0.5 0 0 0;
     0 1.5 0.5 0 0 0;
     0 1 3 0 0 0];
  % matrix for L8 (see part a)

x_dot = A * x;

q_ddot = x_dot(4:6);