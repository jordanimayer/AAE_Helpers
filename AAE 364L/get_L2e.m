%%%%%
% Jordan Mayer
% AAE 364L
% Lab 06
%
% get_L2e:
%   Calculate L2[e] performance criterion for some controller.
%
% Inputs:
%   data: name of data structure as string (e.g. 'errtheta_pitch_benny')
%
% Outputs:
%   L2e: that performance criterion for this data
%%%%%

function [L2e] = get_L2e(data)
    t = data.time;  % sec
    t_f = t(end);  % sec
    theta = data.signals.values(:,1);  % deg
    theta_ref = data.signals.values(:,2);  % deg
    theta_err = theta - theta_ref;  % deg
    psi = errtheta_yaw_benny.signals.values(:,1);  % deg
    psi_ref = errtheta_yaw_benny.signals.values(:,2);  % deg
    psi_err = psi - psi_ref;  % deg
    
    L2e = sqrt((1/tf)*trapz(t,abs(e).^2));  % L2[e] error
end