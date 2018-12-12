%%%%%
% Jordan Mayer
% AAE 532
%
% get_sma_lambert_elliptic:
%   Iterate to get semimajor axis using lambert's equation.
%
% Inputs:
%   type: transfer type, must be one of: 1A, 2A, 1B, 2B
%   mu: gravitational parameter of orbited body, km^3/s^2
%   TOF: time of flight, s
%   c: chord length, km
%   s: semiperimeter, km
%   a_range: range of semimajor axis values to iterate over, km
%            [a_min, a_max]
%
% Outputs:
%   a: semimajor axis, km
%   alpha: alpha from Lambert's equation in proper quadrant, rad
%   beta: beta from Lambert's equation in proper quadrant, rad

function [a, alpha, beta] = get_sma_lambert_elliptic(...
    type, mu, TOF, c, s, a_range)

    alpha = 123; beta = 123;  % rad

    for a = linspace(a_range(1), a_range(2), 10^6)
        
       alpha_0 = 2*asin(sqrt(s/(2*a)));  % rad
       beta_0 = 2*asin(sqrt((s-c)/(2*a)));  % rad
       
       switch type
           case '1A'
               alpha = alpha_0;
               beta = beta_0;
           case '2B'
               alpha = 2*pi - alpha_0;
               beta = -beta_0;
           case '1B'
               alpha = 2*pi - alpha_0;
               beta = beta_0;
           case '2A'
               alpha = alpha_0;
               beta = -beta_0;
           otherwise
               msg = ['Invalid type! Must be one of: ''1A'', ''2A'',' ...
                      '''1B'', ''2A'''];
               error(msg);
       end
       
       TOF_guess = sqrt(a^3 / mu) * ...
                   ((alpha - sin(alpha)) - (beta - sin(beta)));  % s
       err = abs(TOF - TOF_guess);
       
       if err < TOF * 10^(-6)
           break;
       end
    end
    
    if a == a_range(2)
        msg = 'Could not find semimajor axis within range!';
        TOF_mini = TOF_guess
        error(msg);
    end
end

