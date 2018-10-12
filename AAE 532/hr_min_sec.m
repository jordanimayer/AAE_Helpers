%%%%%
% Jordan Mayer
%
% day_hr_min_sec:
%   Convert time in seconds to time in hours, minutes, seconds for
%   more intuitive results.
%
%  Inputs:
%    time: time, s
%
%  Outputs:
%    hr: hours component
%    min: minutes component
%    sec: seconds component
%%%%%

function [hr, min, sec] = hr_min_sec(time)

hr = floor(time/(60*60));
min = floor((time - hr*60*60)/60);
sec = round(mod(time,60));