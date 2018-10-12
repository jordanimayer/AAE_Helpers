%%%%%
% Jordan Mayer
%
% day_hr_min_sec:
%   Convert time in seconds to time in days, hours, minutes, seconds for
%   more intuitive results.
%   Assumes 24 hour day.
%
%  Inputs:
%    time: time, s
%
%  Outputs:
%    day: days component0
%    hr: hours component
%    min: minutes component
%    sec: seconds component
%%%%%

function [day, hr, min, sec] = day_hr_min_sec(time)

day = floor(time/(24*60*60));
hr = floor((time - day*24*60*60)/(60*60));
min = floor((time - day*24*60*60 - hr*60*60)/60);
sec = round(mod(time,60));