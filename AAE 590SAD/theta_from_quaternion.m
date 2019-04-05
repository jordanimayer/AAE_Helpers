%%%%%
% Jordan Maye
% AAE 590: Spacecraft Attitude Dynamics
%
% theta_from_quaternion:
%   Calculate rotation angle using quaternion from initial to final frame.
%
% Inputs:
%   q: quaternion from initial to final frame, vector componentes expressed
%      in either frame
%   theta_max_diff: maximum difference for two angles to be considered
%                   "equal", deg
%
% Outputs:
%   theta: rotation angle for simple rotation
%%%%%

function theta = theta_from_quaternion(q, theta_max_diff)
    norm_q_vec = norm(q(1:3));
    theta_c1 = 2*acosd(q(4));  % cosine option one, deg
    theta_c2 = -theta_c1;  % cosine option two, deg
    theta_s1 = 2*asind(norm_q_vec);  % sine option one, deg
    theta_s2 = 180 - theta_s1;  % sine option two, deg

    theta_c1 = bound_angle_180(theta_c1);
    theta_c2 = bound_angle_180(theta_c2);
    theta_s1 = bound_angle_180(theta_s1);
    theta_s2 = bound_angle_180(theta_s2);

    if ((abs(theta_c1 - theta_s1) < theta_max_diff) || ...
        (abs(theta_c1 - theta_s2) < theta_max_diff))
        theta = theta_c1;
    elseif ((abs(theta_c2 - theta_s1) < theta_max_diff) || ...
            (abs(theta_c2 - theta_s2) < theta_max_diff))
        theta = theta_c2;
    else
        error('Quadrant check for theta failed! No matching angles.');
    end
end