%%%%%
% Jordan Mayer
% AAE 532
%
% Saturn_flyby:
%   Calculate outbound velocity for Saturn flyby given arrival conditions
%   and closest approach distance.
%
% Inputs:
%   r_p: closest approach distance, km
%   v_infSaturn_minus: inbound velocity vector relative to Saturn, km/s,
%                      x-y-z coords
%   v_Saturn: Saturn velocity vector, km/s, x-y-z coords
%   del: turn angle, rad
%
% Outputs:
%   v_infSaturn_plus: outbound velocity vector relative to Saturn, km/s,
%                     x-y-z coords
%                     both options, as array: [option1, option2]
%%%%%

function [v_infSaturn_plus] = Saturn_flyby(r_p, v_infSaturn_minus, ...
    v_Saturn, del)

    x_minus = v_infSaturn_minus(1)  % km/s
    y_minus = v_infSaturn_minus(2)  % km/s
    z_minus = v_infSaturn_minus(3)  % km/s
    v_infSaturn = norm(v_infSaturn_minus)  % km/s
    v_dot = v_infSaturn^2 * cos(del)  % km/s
    
    fprintf('\n');
    
    a = (x_minus/y_minus)^2 + 1
    b = 2*(v_dot/y_minus)*(-x_minus/y_minus)
    c = (v_dot/y_minus)^2 - v_infSaturn^2
    
    fprintf('\nOption 1:\n\n');
    
    x_plus_1 = (-b + sqrt(b^2 - 4*a*c))/(2*a)  % km/s
    y_plus_1 = v_dot/y_minus - (x_minus/y_minus)*x_plus_1  % km/s
    v_infSaturn_plus_1 = [x_plus_1; y_plus_1; 0]  % km/s, x-y-z coords
    
    err_mag_1 = abs(norm(v_infSaturn_plus_1) - v_infSaturn);
    err_dot_1 = abs(dot(v_infSaturn_minus, v_infSaturn_plus_1) - v_dot);
    if (err_mag_1 > 0.0001 || err_dot_1 > 0.0001)
        error('Calculated v_infSaturn_plus_1 doesn''t check out!');
    end
    
    fprintf('\nOption 2:\n\n');
    
    x_plus_2 = (-b - sqrt(b^2 - 4*a*c))/(2*a)  % km/s
    y_plus_2 = v_dot/y_minus - (x_minus/y_minus)*x_plus_2  % km/s
    v_infSaturn_plus_2 = [x_plus_2; y_plus_2; 0]  % km/s, x-y-z coords
    
    err_mag_2 = abs(norm(v_infSaturn_plus_2) - v_infSaturn);
    err_dot_2 = abs(dot(v_infSaturn_minus, v_infSaturn_plus_2) - v_dot);
    if (err_mag_2 > 0.0001 || err_dot_2 > 0.0001)
        error('Calculated v_infSaturn_plus_1 doesn''t check out!');
    end
    
    fprintf('\n');
    
    v_infSaturn_plus = [v_infSaturn_plus_1, v_infSaturn_plus_2];
        % both options
end
