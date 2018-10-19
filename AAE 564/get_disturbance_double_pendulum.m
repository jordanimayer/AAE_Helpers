%%%%%
% Jordan Mayer
% AAE 564
% HW 08
%
% get_disturbance_double_pendulum:
%   Calculate disturbance, w, for double pendulum system
%
%   Inputs:
%     alpha
%     omega
%     t (time, s)
%
%   Output:
%     w
%%%%%

function [w] = get_disturbance_double_pendulum(alpha, omega, t)
  w = exp(alpha*t) * sin(omega*t);
end