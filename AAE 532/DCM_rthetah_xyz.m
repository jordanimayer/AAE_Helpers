%%%%%
% Jordan Mayer
% AAE 532
%
% DCM_rthetah_xyz:
%   Calculate DCM from r-theta-h coordinates to x-y-z coordinates
%   (using equatorial plane).
%
%   inputs:
%       raan:  right ascension of ascending node (deg)
%       i: inclination (deg)
%       omega: argument of periapsis (deg)
%       thetaStar: true anomaly (deg)
%
%   output:
%       DCM: transformation matrix from r-theta-h coordinates to x-y-z
%            coordinates (using equatorial plane)
%%%%%

function [DCM] = DCM_rthetah_xyz(raan, i, omega, thetaStar)
    % get theta, keeping between -180 and 180 degrees
    theta = omega + thetaStar;
    if (theta > 180)
        theta = theta - 360;
    end
    if (theta < -180)
        theta = theta + 360;
    end
    
    % calculate terms for DCM elements
    sRaan = sind(raan);
    cRaan = cosd(raan);
    sTheta = sind(theta);
    cTheta = cosd(theta);
    sI = sind(i);
    cI = cosd(i);
    
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
    if  abs(norm(r_hat) - 1) > 0.00001 ||...
        abs(norm(theta_hat) - 1) > 0.00001 ||...
        abs(norm(h_hat) - 1) > 0.00001 ||...
        abs(norm(x_hat) - 1) > 0.00001 ||...
        abs(norm(y_hat) - 1) > 0.00001 ||...
        abs(norm(z_hat) - 1) > 0.00001
        
        error('Not all rows/columns are unit vectors!');
    end
        
end

