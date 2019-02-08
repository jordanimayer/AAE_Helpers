%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% kron_delta:
%   Calculate Kronecker delta for two numbers.
%
% Inputs:
%   i, j: numbers
%
% Outputs:
%   delta: Kronecker delta of i and j
%%%%%

function [delta] = kron_delta(i, j)
    delta = (i == j);
end