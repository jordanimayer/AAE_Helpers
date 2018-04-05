%%%%%
% AAE 364
% HW 07
% Jordan Mayer
%
% Given a control system, find the poles such that the dominant poles have
% a damping ratio of 0.5.
% This is accomplished by solving iteratively for the positive
% complex pole using the angle condition.
%%%%%

%%% Problem 1: B-6-6 %%%

zeta = 0.5;                     % damping ratio (known)
wn = linspace(0, 10, 10^6);     % possible values of natural frequency
                                % (unknown)
                                
% Iterate over natural frequency values
% For each natural frequency, calculate the corresponding
% positive complex pole, s.
% Then calculate L(s) and the angle of L(s) and compare to -180 degrees.
% Once an s has been found such that L(s) is equal to -180 degrees
% (+/- 0.001), stop iterating and display the positive complex pole
% and the corresponding angle of L(s).

for k=1:10^6
    s = -zeta*wn(k) + 1i*wn(k)*sqrt(1-zeta^2);
    L = (s+9)/(s*(s^2+4*s+11));
    if (abs(rad2deg(angle(L)) + 180) < 0.001)
        break;
    end
end
disp(s);                    % output: -1.5000 + 2.5981i
disp(rad2deg(angle(L)));    % output: -179.9994

% Note that, due to symmetry about the Real axis, the negative complex
% pole must be -1.5000 - 2.5981i.
% The real pole is solved for in part c) after solving for the
% corresponding k using the characteristic equation.