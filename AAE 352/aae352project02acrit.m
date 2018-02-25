% AAE 352 - Project 02
% Solve for a_crit
% Jordan Mayer

K_IC = 26;  %K_IC, MPa*m^(1/2)
w = 0.023;  %w, m
sigmaRef = 34.69;   %sigmaRef, MPa

a_crit = 0; %a_crit, m

% Find a_crit iteratively
for a=linspace(0,0.5,100000)
    x = a/w;
    beta = 30.795*x^4 - 51.44*x^3 + 29.462*x^2 - 6.2025*x + 2.0791;
    K = beta*sigmaRef*sqrt(pi*a);
    if (K >= K_IC)
        disp(K);
        a_crit = a
        break;
    end
end

C = 10^(-9.1141);   % constant for da/dN calculation
m = 2.4881;         % constant for da/dN calculation

a_prev = 0; % crack size at last inspection, m

% Find a_prev iteratively
for a=linspace(0,0.5,100000)
    x = a/w;
    beta = 30.795*x^4 - 51.44*x^3 + 29.462*x^2 - 6.2025*x + 2.0791;
    K = beta*sigmaRef*sqrt(pi*a);
    dadN = C*K^m;
    acrit = a + dadN*11000;
    if (acrit >= a_crit)
        disp(acrit);
        a_prev = a
        break;
    end
end