%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% theta_from_sin_cos:
%   Calculate angle given cosine and sine.
%
% Inputs:
%   sin_theta: sine of angle
%   cos_theta: cosine of angle
%   theta_max_diff: maximum difference for two angles to be considered
%                   equal (used for quadrant checks)
%
% Outputs:
%   theta: angle, deg
%%%%%

function [theta] = theta_from_sin_cos(sin_theta, cos_theta, theta_max_diff)
    theta_s1 = asind(sin_theta);
    theta_s2 = 180 - theta_s1;
    theta_c1 = acosd(cos_theta);
    theta_c2 = -theta_c1;
    
    theta_s1 = bound_angle_360(theta_s1);
    theta_s2 = bound_angle_360(theta_s2);
    theta_c1 = bound_angle_360(theta_c1);
    theta_c2 = bound_angle_360(theta_c2);
    
    if ((abs(theta_c1 - theta_s1) < theta_max_diff) || ...
        (abs(theta_c1 - theta_s2) < theta_max_diff))
        theta = theta_c1;
    elseif ((abs(theta_c2 - theta_s1) < theta_max_diff) || ...
            (abs(theta_c2 - theta_s2) < theta_max_diff))
        theta = theta_c2;
    else
        disp(sin_theta);
        disp(cos_theta);
        disp(theta_max_diff);
        error('Quadrant check for theta failed! No matching angles.');
    end
end