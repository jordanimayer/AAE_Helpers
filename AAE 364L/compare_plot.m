%%%%%
% Jordan Mayer
% AAE 364L
% Lab 06
%
% compare_plot:
%   Do that.
%
% Inputs:
%   data: data structure with .time and .signals.values (two 'o dem)
%   ylab: y-axis label
%   ttl: title
%%%%%

function [] = compare_plot(data1, data2, ylab, ttl)
    t1 = data1.time;  % sec
    t2 = data2.time;  % sec
    angle1 = data1.signals.values;  % deg
    angle_ref = angle1(:,2);  % deg
    angle1 = angle1(:,1);  % deg
    angle2 = data2.signals.values(:,1);  % deg
    
    figure();
    plot(t1, angle_ref, '-k'); hold on;
    plot(t1, angle1, '-b'); plot(t2, angle2, '-r');
    title(ttl);
    xlabel('Time, sec'); ylabel(ylab); grid on;
    legend('Reference', 'Error Feedback Response', ...
           'State Feedback Response');
end