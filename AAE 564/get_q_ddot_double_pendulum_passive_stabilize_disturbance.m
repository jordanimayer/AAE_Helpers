%%%%%
% Jordan Mayer
% AAE 564
% HW 08
%
% get_q_ddot_double_pendulum_passive_stabilize_disturbance:
%   Get q_ddot for double-pendulum problem, using equation M*q_ddot = F
%
%   Inputs:
%     m: [m0, m1, m2]'
%     l: [l1, l2]'
%     u
%     q: [y, theta1, theta2]'
%     q_dot: [y_dot, theta1_dot, theta2_dot]'
%     k
%     c
%     alpha
%     omega
%     t
%
%   Output:
%     q_ddot: [y_ddot, theta1_ddot, theta2_ddot]
%%%%%

function q_ddot = ...
    get_q_ddot_double_pendulum_passive_stabilize_disturbance(...
      m, len, q, q_dot, k, c, alpha, omega, t)
g = 1;  % constant defined in problem
% extract inputs for easier readability
m0 = m(1);
m1 = m(2);
m2 = m(3);
len1 = len(1);
len2 = len(2);
y = q(1);
y_dot = q_dot(1);
theta1 = q(2);
theta1_dot = q_dot(2);
theta2 = q(3);
theta2_dot = q_dot(3);
a = exp(alpha*t);
w = sin(omega*t);
u = -k*y - c*y_dot + w;

% define intermittent terms to populate matrices
mlcosth1 = m1*len1*cos(theta1);
mlcosth2 = m2*len2*cos(theta2);
mlsinth1 = m1*len1*sin(theta1);
mlsinth2 = m2*len2*sin(theta2);

% define matrices
M = [sum(m), -mlcosth1, -mlcosth2;
     -mlcosth1, m1*len1*len1, 0;
     -mlcosth2, 0, m2*len2*len2];
F = [u - mlsinth1*theta1_dot*theta1_dot - mlsinth2*theta2_dot*theta2_dot;
     -g*mlsinth1;
     -g*mlsinth2];

q_ddot = M\F;