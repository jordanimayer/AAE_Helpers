%%%%%
% Jordan Mayer
% AAE 564
%
% get_A_B_double_pendulum:
%   Obtain matrices A and B for state-space representation of linearized
%   double-pendulum system.
%
% Inputs:
%   P: [m0 m1 m2 len1 len2 g]
%   E: [y_e theta1_e theta2_e] (thetas in deg)
%
% Output:
%   A and B matrices
%%%%%

function [A, B] = get_A_B_double_pendulum(P,E)
    m0 = P(1); m1 = P(2); m2 = P(3); len1 = P(4); len2 = P(5); g = P(6);
    y_e = E(1); theta1_e = E(2); theta2_e = E(3);
    
    M = zeros(3,3); K = zeros(3,3); phi = [1 0 0]';
    
    M(1,1) = m0 + m1 + m2;
    M(1,2) = -m1*len1*cosd(theta1_e);
    M(1,3) = -m2*len2*cosd(theta2_e);
    M(2,1) = M(1,2);
    M(2,2) = m1*len1^2;
    M(3,1) = M(1,3);
    M(3,3) = m2*len2^2;
    
    K(2,2) = -m1*len1*g*cosd(theta1_e);
    K(3,3) = -m2*len2*g*cosd(theta2_e);
    
    A = [zeros(3,3) eye(3,3); inv(M)*K zeros(3,3)];
    B = [0; 0; 0; inv(M)*phi];
end

