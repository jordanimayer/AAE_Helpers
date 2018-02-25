% AAE 334 - HW 08
% Jordan Mayer

% Problem 3

M = linspace(0.0, 2.0, 100);
prat = zeros(3,100);    % pressure ratio (po2/p1)
gamma = [1.67, 1.4, 1.044];

% Calculate pressure ratios
for i=1:3
    for j=1:100
        if (M(j) < 1) %isentropic
            prat(i,j) = (1+(gamma(i)-1)/2 * M(j)^2)^(gamma(i)/(gamma(i)-1));
        else %supersonic
            prat(i,j) = (((gamma(i)+1)^2 * M(j)^2)/(4*gamma(i)*M(j)^2 - 2*(gamma(i)-1)))^(gamma(i)/(gamma(i)-1));
            prat(i,j) = prat(i,j)*((1-gamma(i)+2*gamma(i)*M(j)^2)/(gamma(i)+1));
        end
    end
end

% Create and format plots
figure(1)
plot(M, prat(1, 1:100), '-r'); hold on;
plot(M, prat(2, 1:100), '-b'); hold on;
plot(M, prat(3, 1:100), '-k'); grid on;
xlim([0, 2]);ylim([0, 7]);
title('Pitot Tube Pressure Prediction');
xlabel('M');ylabel('p_{o2}/p_{1}');
legend('Argon', 'Air', 'Octane');


% Problem 4
To = 300; % stagnation temperature, K
Ae = 0.0001; % exit area, m^2
po = 100000; % stagnation pressure, Pa
pe = linspace(0, po, 1000); % exit pressure, Pa
mdot = zeros(1, 1000); % mass flow rate, kg/s
mdotnd = zeros(1, 1000);    % nondimensional mass flow rate
pend = zeros(1, 1000);  % nondimensional exit pressure
gamma = 1.4; % specific heat ratio
R = 287; % ideal gas constant for air
ao = sqrt(gamma*R*To);  % stagnation sound speed, m/s
pe3 = po*(2/(gamma+1))^(gamma/(gamma-1));   % sonic pressure

% Calculate mass flow rates
for j=1:1000
    if (pe(j) >= pe3)
        Me = (2/(gamma-1)*((po/pe(j))^((gamma-1)/gamma)-1))^(1/2); % exit Mach
        rhoo = po/(R*To);   % stagnation density, kg/m^3
        rhoe = rhoo*(1+(gamma-1)/2*Me^2)^(-1/(gamma-1)); % exit density, kg/m^3
        Te = To*(1+(gamma-1)/2*Me^2)^(-1); % exit temperature, K
        ue = Me*sqrt(gamma*R*Te); % exit velocity, m/s
        mdot(j) = rhoe*ue*Ae; % UNCHOKED mass flow rate, kg/s
    else
        pstar = po*(2/(gamma+1))^(gamma/(gamma-1)); % sonic (choked) pressure
        Tstar = To*(2/(gamma+1)); % sonic (choked) temperature
        mdot(j) = pstar*Ae*sqrt(gamma/(R*Tstar));   % CHOKED mass flow rate, kg/s
    end

    mdotnd(j) = mdot(j)/(rhoo*Ae*ao);
    pend(j) = pe(j)/po;
end

% Plot results
figure(2)
plot(pe, mdot, '-k'); grid on;
title('Converging Nozzle Flow');
xlabel('Exit pressure, Pa'); ylabel('Mass flow rate, kg/s');
figure(3)
plot(pend, mdotnd, '-k'); grid on;
title('Convergin Nozzle Flow (Nondimensional)');
xlabel('Nondimensional exit pressure'); 
ylabel('Nondimensional mass flow rate');