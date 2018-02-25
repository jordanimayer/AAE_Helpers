% AAE 334 - HW 11
%
% Jordan Mayer

% Problem 1

gamma = 1.4;    % specific heat ratio for air
Mcr = linspace(0, 1, 1000);     % critical Mach number to find intersection
Cp0 = linspace(-15, -0.001, 1000);   % stagnation pressure coefficient
Mcr_PG = zeros(1,1000);     % critical Mach numbers from Prandtl-Glauert Rule
Mcr_KT = zeros(1,1000);     % critical Mach numbers from Karman-Tsien Rule
Mcr_L  = zeros(1,1000);     % critical Mach numbers from Laitone Rule

for i=1:1000    % Cp0 values (for plot)
    for j=1:1000    % Mcr values (for finding intersection)
        % Calculate Cpcr from flow relation
        Cpcr = (1 + (gamma-1)/2 * Mcr(j)^2)/((gamma+1)/2);
        Cpcr = Cpcr^(gamma/(gamma-1)) - 1;
        Cpcr = Cpcr * (2/(gamma*Mcr(j)^2));

        % Calculate Cpcr from Prandtl-Glauert Rule
        Cpcr_PG = Cp0(i)/(sqrt(1-Mcr(j)^2));
        
        % Calculate Cpcr from Karman-Tsien Rule
        Cpcr_KT = Cp0(i)/(sqrt(1-Mcr(j)^2) + (Mcr(j)^2 / (1+sqrt(1-Mcr(j)^2))) * Cp0(i)/2);
        
        % Calculate Cpcr from Laitone Rule
        Cpcr_L = Cp0(i)/(sqrt(1-Mcr(j)^2) + (Mcr(j)^2 * (1+(gamma-1)/2 * Mcr(j)^2) / (2*sqrt(1-Mcr(j)^2))) * Cp0(i));
    
        if (Cpcr >= Cpcr_PG && Mcr_PG(i) == 0) 
            % Cpcr has just intersected Cpcr_PG
            Mcr_PG(i) = Mcr(j);
        end
        if (Cpcr >= Cpcr_KT && Mcr_KT(i) == 0)
            % Cpcr has just intersected Cpcr_KT
            Mcr_KT(i) = Mcr(j);
        end
        if (Cpcr >= Cpcr_L && Mcr_L(i) == 0)
            % Cpcr has just intersected Cpcr_L
            Mcr_L(i) = Mcr(j);
        end
        
        if (Mcr_PG(i) * Mcr_KT(i) * Mcr_L(i) ~= 0)
            % Cpcr has intersected with all other Cpcrs
            % --> go to next Cp0 value
            break;
        end
    end
end

% Create plot
plot(Cp0, Mcr_PG); hold on;
plot(Cp0, Mcr_KT);
plot(Cp0, Mcr_L);

% Format plot
legend('Prandtl-Glauert Rule', 'Karman-Tsien Rule', 'Laitone Rule');
title('Prediction of Critical Mach Number');
xlabel('C_{p0}'); ylabel('M_{cr}');