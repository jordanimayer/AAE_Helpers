%%%%%
% Jordan Mayer
% AAE 590: Spacecraft Attitude Dynamics
%
% get_K_from_I:
%   Calculate K variables from principal moments of inertia.
%
% Inputs:
%   I: Principal moments of inerta, [I_1, I_2, I_3], consistent units
%
% Outputs:
%   K: K variables, [K_1, K_2, K_3]
%%%%%

function K = get_K_from_I(I)
    K = zeros(1,3);
    
    K(1) = (I(2)-I(3))/I(1);
    K(2) = (I(3)-I(1))/I(2);
    K(3) = (I(1)-I(2))/I(3);
end