%%%%%
% Jordan Mayer
% AAE 532
%
% rthetah_to_xyz:
%   Use DCM to convert from r-theta-h coordinates to x-y-z coordinates
%   (using equatorial plane).
%
%   inputs:
%       v_rth: vector in r-theta-h coordinates
%       i: inclination (rad)
%       raan: right ascension of ascending node (rad)
%       theta: true anomaly plus argument of periapsis (rad)
%
%   output:
%       v_xyz: same vector in x-y-z coordinates
%%%%%

function [v_xyz] = rthetah_to_xyz(v_rth, i, raan, theta)
    
    % calculate terms for DCM elements
    sI = sin(i);
    cI = cos(i);
    sRaan = sin(raan);
    cRaan = cos(raan);
    sTheta = sin(theta);
    cTheta = cos(theta);
    
    % calculate DCM element-by-element
    DCM = zeros(3,3);
    DCM(1,1) = cRaan*cTheta - sRaan*cI*sTheta;
    DCM(1,2) = -cRaan*sTheta - sRaan*cI*cTheta;
    DCM(1,3) = sRaan*sI;
    DCM(2,1) = sRaan*cTheta + cRaan*cI*sTheta;
    DCM(2,2) = -sRaan*sTheta + cRaan*cI*cTheta;
    DCM(2,3) = -cRaan*sI;
    DCM(3,1) = sI*sTheta;
    DCM(3,2) = sI*cTheta;
    DCM(3,3) = cI;
    
    % check DCM rows/columns
    r_hat = DCM(:,1);
    theta_hat = DCM(:,2);
    h_hat = DCM(:,3);
    x_hat = DCM(1,:)';
    y_hat = DCM(2,:)';
    z_hat = DCM(3,:)';
    if  abs((norm(r_hat) * norm(theta_hat) * norm(h_hat) * ...
             norm(x_hat) * norm(y_hat) * norm(z_hat)) - 1) > 0.000001
        
        error('Not all rows/columns are unit vectors!');
    end
    
    % convert vector
    v_xyz = DCM * v_rth;
end
