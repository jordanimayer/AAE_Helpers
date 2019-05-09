%%%%%
% Jordan Mayer
% AAE 364L
% Lab 06
%
% create_plot:
%   Do that.
%
% Inputs:
%   data: data structure with .time and .signals.values
%   ylab: x-axis label
%   ttl: title
%%%%%

function [] = create_plot(data, ylab, ttl)
    t = data.time;  % sec
    angle = data.signals.values;  % deg
    angle_ref = angle(:,2);  % deg
    angle = angle(:,1);  % deg
    
    figure();
    plot(t, angle_ref, '-b'); hold on; plot(t, angle, '-r');
    title(ttl);
    xlabel('Time, sec'); ylabel(ylab); grid on;
    legend('Reference', 'Actual');
end