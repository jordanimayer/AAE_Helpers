%%%%%
% Jordan Mayer
% AAE 532
%
% zero_to_twoPi:
%   Ensure angle is expressed as a number between 0 and 2pi radians.
%
% Inputs:
%   theta: angle, rad
%
% Outputs:
%   theta_corrected: angle, rad, between 0 and 2pi.
%%%%%

function [theta_corrected] = zero_to_twoPi(theta)
    theta_corrected = theta;
    
    while (theta_corrected < 0)
        theta_corrected = theta_corrected + 2*pi;
    end
    
    while (theta_corrected > 2*pi)
        theta_corrected = theta_corrected - 2*pi;
    end
end
