%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% bound_angle_180:
%   Transform angle to be between -180 and 180 degrees.
%
% Inputs:
%   theta: angle, in deg
%
% Outputs:
%   theta_bound: angle, transformed by adding/subtracting 360 deg to be
%                  between -180 and 180 degrees
%%%%%

function theta_bound = bound_angle_180(theta)
    while (theta > 180)
        theta = theta - 360;
    end
    while (theta < -180)
        theta = theta + 360;
    end
    
    theta_bound = theta;
end