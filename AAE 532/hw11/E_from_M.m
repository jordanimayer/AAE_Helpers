%%%%%
% Jordan Mayer
% AAE 532
%
% E_from_M:
%   Iterate to calculate eccentric anomaly from mean anomaly for elliptical
%   orbit.
%
% Inputs:
%   e: eccentricity of orbit
%   M: mean anomaly, rad
%
% Outputs:
%   E: eccentric anomaly, rad

function [E] = E_from_M(e, M)
    E = 4321;  % indicates not found
    
    for E_guess = linspace(0, 2*pi, 1000000)
        M_result = E_guess - e*sin(E_guess);
        if abs(M - M_result) < 0.00005
            E = E_guess;
            break;
        end
    end
end
