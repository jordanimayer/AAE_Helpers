%%%%%
% Jordan Mayer
% AAE 364L
% Lab 06
%
% plot_volt:
%   Do that.
%
% Inputs:
%   data: data structure with .time and .signals.values
%   ttl: title
%%%%%

function [] = plot_volt(data, ttl)
    t = data.time;  % sec
    v = data.signals.values;  % V
    vp = v(:,1);  % V
    vy = v(:,2);  % V
    
    figure();
    plot(t, vp, '-b'); hold on; plot(t, vy, '-r');
    title(ttl);
    xlabel('Time, sec'); ylabel('Input Voltage, V'); grid on;
    legend('Pitch', 'Yaw');
end