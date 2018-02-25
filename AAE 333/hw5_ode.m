function [ dvdt ] = hw5_ode( t, v )
% AAE 333: Fluid Mechanics - HW 5, Problem 3
% This function calculates dv/dt for a small rocket

dvdt = (1957/(-5.76*t)) - 0.00146 * ((v^2)/(-5.76*t)) - 9.81;

end

