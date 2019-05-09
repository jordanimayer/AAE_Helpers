%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% theta_from_cos_sin_unbound:
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

function [theta] = theta_from_sin_cos_unbound(sin_theta, cos_theta, ...
                                              theta_max_diff)
    theta_s1 = asind(cos_theta);
    theta_s2 = 180 - theta_s1;
    theta_c1 = acosd(sin_theta);
    theta_c2 = -theta_c1;
    
    theta_s1_b = bound_angle_360(theta_s1);
    theta_s2_b = bound_angle_360(theta_s2);
    theta_c1_b = bound_angle_360(theta_c1);
    theta_c2_b = bound_angle_360(theta_c2);
    
    if ((abs(theta_c1_b - theta_s1_b) < theta_max_diff) || ...
        (abs(theta_c1_b - theta_s2_b) < theta_max_diff))
        theta = theta_c1;
    elseif ((abs(theta_c2_b - theta_s1_b) < theta_max_diff) || ...
            (abs(theta_c2_b - theta_s2_b) < theta_max_diff))
        theta = theta_c2;
    else
        disp(sin_theta);
        disp(cos_theta);
        disp(theta_max_diff);
        error('Quadrant check for theta failed! No matching angles.');
    end
end