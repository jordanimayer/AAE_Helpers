% AAE 352 - HW 10
%
% Jordan Mayer

% Problem 1

syms a;
m = 3.2214;
C = 10^(-11.06);
sigmaMax = 125; % NOTE: keep units in MPa, because m and C were found with
                % deltaK in MPa*sqrt(m)
a_i = 0.01*10^(-3); % NOTE: convert units to m, because m and C were found
                    % with deltaK in MPa*sqrt(m)
a_f = 1.00*10^(-3);
beta = 1.12;

N = int(1/(C*(beta*sigmaMax*sqrt(pi*a))^m), a, a_i, a_f);
N = double(N)