%%%%%
% Jordan Mayer
% AAE 364L
% Lab 06
%
% get_errs:
%   Calculate L2[e] and eMss performance criterion for some controller.
%
% Inputs:
%   data: data structure with .time and .signals.values
%
% Outputs:
%   L2e: that performance criterion for this data
%   eMss: that one
%%%%%

function [L2e, eMss] = get_errs(data)
    t = data.time;  % sec
    t_f = t(end);  % sec
    angle = data.signals.values(:,1);  % deg
    angle_ref = data.signals.values(:,2);  % deg
    e = angle - angle_ref;  % deg
    
    L2e = sqrt((1/t_f)*trapz(t,abs(e).^2));  % L2[e] error
    
    [ind, tmp] = find(t >= 30,1);
    if t(ind) ~= 30
        error('t(ind) is not 30!');
    end
    
    eMss = max(abs(e(ind:end)));
end