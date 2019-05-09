%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% bound_angle_360:
%   Transform angle to be between 0 and 360 degrees.
%
% Inputs:
%   theta: angle, in deg
%
% Outputs:
%   theta_bound: angle, transformed by adding/subtracting 360 deg to be
%                  between 0 and 360 degrees
%%%%%

function theta_bound = bound_angle_360(theta)
    while (theta > 360)
        theta = theta - 360;
    end
    while (theta < 0)
        theta = theta + 360;
    end
    
    theta_bound = theta;
end