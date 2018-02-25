%%%
% AAE 339 - HW 03
% Jordan Mayer
%%%

%%% Problem 1 %%%

% Given a set of expansion ratio values, plot:
% static temperature, pressure, and density
% velocity, sound speed, and Mach number

% Initialize variables
    % Known
expansion = [2 5 10 20];    % expansion ratios
T0 = 3000;                  % stagnation temperature, K
p0 = 10.0*10^6;             % stagnation pressure, Pa
gamma = 1.2;                % specific heat ratio, dimensionless
M = [2.055 2.785 3.278 3.759];  % Mach number, from tables
R = 287.1;                  % ideal gas constant (same as air), J/(kg*K)

    % To calculate
T = zeros(1,4);             % static temperature, K
p = zeros(1,4);             % static pressure, Pa
rho = zeros(1,4);           % static density, kg/m^3
u = zeros(1,4);             % velocity, m/s
a = zeros(1,4);             % sound speed, m/s

% Calculate values
for k = 1:4   % for each expansion ratio
    % solve for temperature
    T(k) = T0/(1+(1/2)*(gamma-1)*M(k)^2);
    p(k) = p0*((T(k)/T0)^(gamma/(gamma-1)));
    rho(k) = p(k)/(R*T(k));
    a(k) = sqrt(gamma*R*T(k));
    u(k) = M(k)*a(k);
end

% Create and format plots
% close all;
% 
% figure(1);plot(expansion, T);hold on;plot(expansion, T, 'ro');
% title('Rocket Nozzle Analysis: Temperature');
% xlabel('A/A*');ylabel('T, K');grid on;
% 
% figure(2);plot(expansion, p*10^(-6));hold on;plot(expansion, p*10^(-6), 'ro');
% title('Rocket Nozzle Analysis: Pressure');
% xlabel('A/A*');ylabel('p, MPa');grid on;
% 
% figure(3);plot(expansion, rho);hold on;plot(expansion, rho, 'ro');
% title('Rocket Nozzle Analysis: Density');
% xlabel('A/A*');ylabel('rho, kg/m^3');grid on;
% 
% figure(4);plot(expansion, u);hold on;plot(expansion, u, 'ro');
% title('Rocket Nozzle Analysis: Velocity');
% xlabel('A/A*');ylabel('u, m/s');grid on;
% 
% figure(5);plot(expansion, a);hold on;plot(expansion, a, 'ro');
% title('Rocket Nozzle Analysis: Sound Speed');
% xlabel('A/A*');ylabel('a, m/s');grid on;
% 
% figure(6);plot(expansion, M);hold on;plot(expansion, M, 'ro');
% title('Rocket Nozzle Analysis: Mach Number');
% xlabel('A/A*');ylabel('M');grid on;


%%% Problem 2 %%%

% Initialize variables
    % Known
        % Air properties
gamma = 1.4;            % specific heat ratio
R = 287.1;              % ideal gas constant, J/(kg*K)
cp = 1.0*10^3;          % constant-pressure specific heat capacity
    % State 1
T01 = 560.7;            % stagnation temperature, K
p01 = 7.255*10^5;       % stagnation pressure, Pa
M1 = [0.1 0.3 0.5];     % Mach number

T01rat = [0.04678 0.3469 0.6914];     % temperature ratio, T01/T0*
    
    % To calculate
        % NOTE: ratios obtained using Compressible Aerodynamics Calculator
        % after obtaining values for q, T02, and Trat, 
        % following method explained in handwritten portion.
T02 = zeros(1,4);    % stagnation temperature, K
T02rat = zeros(1,4); % stagnation temperature ratio
p02 = zeros(1,4);    % stagnation pressure, Pa
p2 = zeros(1,4);     % static pressure, Pa
T2 = zeros(1,4);     % static temperature, K
a2 = zeros(1,4);     % sound speed, m/s
u2 = zeros(1,4);     % velocity, m/s
qmax = zeros(1,3);      % specific heat transfer to choke, J/kg

% Calculate unknowns
    % Heat transfer to choke
for k=1:3
    qmax(k) = cp*T01*(1/T01rat(k)-1);
end
    % State 2 if M1 = 0.3
q = linspace(0, qmax(2), 4);     % specific heat transfer, J/kg
for k=1:4
    T02(k) = q(k)/cp + T01;
    T02rat(k) = T02(k)/T01*T01rat(2);
end
fprintf('T02 rat = ');
disp(T02rat);       % display T02/T0* ratios so we can use the
                    % Compressible Aerodynamics Calculator to find
                    % other ratios
    % State 1 properties and ratios, found using CAC for M1 = 0.3
T1 = T01/(1+(1/2)*(gamma-1)*M1(2)^2);  % static temperature, K
p1 = p01*(T1/T01)^(gamma/(gamma-1));    % static pressure, Pa
T01rat = T01rat(2); % stagnation temperature ratio
p01rat = 1.199;     % stagnation pressure ratio
p1rat = 2.131;      % static pressure ratio
T1rat = 0.4089;     % static temperature ratio
    % Mach number and other ratios, found using CAC
M2 = [0.3000 0.4205 0.5680 1];     % Mach number
p02rat = [1.199 1.148 1.087 1];     % stagnation pressure ratio
p2rat = [2.131 1.924 1.653 1];      % static pressure ratio
T2rat = [0.4089 0.6544 0.8819 1];  % static temperature ratio
for k=1:4
    p02(k) = p02rat(k)/p01rat*p01;
    p2(k) = p2rat(k)/p1rat*p1;
    T2(k) = T2rat(k)/T1rat*T1;
    a2(k) = sqrt(gamma*R*T2(k));
    u2(k) = M2(k)*a2(k);
end

% Create plots
% close all;
% x = 'q, MJ/kg';
% q = q*10^(-6);     % conver to MJ/kg
% 
% figure(1);plot(q, p02*10^(-6));hold on;plot(q, p02*10^(-6), 'ro');
% title('Ramjet Analysis: Combustor Stagnation Pressure');
% xlabel(x);ylabel('p_0, MPa');grid on;
% 
% figure(2);plot(q, p2*10^(-6));hold on;plot(q, p2*10^(-6), 'ro');
% title('Ramjet Analysis: Combustor Static Pressure');
% xlabel(x);ylabel('p, MPa');grid on;
% 
% figure(3);plot(q, T2);hold on;plot(q, T2, 'ro');
% title('Ramjet Analysis: Combustor Static Temperature');
% xlabel(x);ylabel('T, K');grid on;
% 
% figure(4);plot(q, a2);hold on;plot(q, a2, 'ro');
% title('Ramjet Analysis: Combustor Sound Speed');
% xlabel(x);ylabel('a, m/s');grid on;
% 
% figure(5);plot(q, u2);hold on;plot(q, u2, 'ro');
% title('Ramjet Analysis: Combustor Exit Velocity');
% xlabel(x);ylabel('u, m/s');grid on;
    
%%% Problem 3 %%%

% Initialize variables
    % Known
        % Oxygen properties
gamma = 1.4;            % specific heat ratio
R = 259.8;              % ideal gas constant, J/(kg*K)
        % State 1 stagnation conditions
T01 = 600;              % stagnation temperature, K
p01 = 20*10^6;          % stagnation pressure, Pa
        % Properties: [<state 1> <state 2> <state 3>]
L = [0 0.05 0.10];              % length, m
M = [0.30 0.3089 0.3186];       % Mach number
p0rat = [2.035 1.983 1.929];    % stagnation pressure ratio, p0/p0*
Trat = [1.179 1.178 1.176];     % temperature ratio, T/T*
prat = [3.619 3.513 3.403];     % pressure ratio, p/p*
urat = [0.3257 0.335 0.3456];   % velocity ratio, u/u*
    % To calculate
p0 = zeros(1,3);        % stagnation pressure, Pa
p = zeros(1,3);         % pressure, Pa
T = zeros(1,3);         % temperature, K
u = zeros(1,3);         % velocity, m/s

% Calculate unknowns
    % State 1
p0(1) = p01; 
T(1) = T01/(1+(1/2)*(gamma-1)*M(1)^2);
p(1) = p01*(T(1)/T01)^(gamma/(gamma-1));
u(1) = M(1)*sqrt(gamma*R*T(1));
    % State 2 and 3
for k=2:3
    p0(k) = p0rat(k)/p0rat(1)*p0(1);
    T(k) = Trat(k)/Trat(1)*T(1);
    p(k) = prat(k)/prat(1)*p(1);
    u(k) = urat(k)/urat(1)*u(1);
end

% Create plots
close all;

figure(1);plot(L*100, p0*10^(-6));hold on;plot(L*100, p0*10^(-6), 'ro');
title('Oxidizer Injection Analysis: Stagnation Pressure');
xlabel('L, cm');ylabel('p_0, MPa');grid on;

figure(2);plot(L*100, p*10^(-6));hold on;plot(L*100, p*10^(-6), 'ro');
title('Oxidizer Injection Analysis: Static Pressure');
xlabel('L, cm');ylabel('p, MPa');grid on;

figure(3);plot(L*100, T);hold on;plot(L*100, T, 'ro');
title('Oxidizer Injection Analysis: Temperature');
xlabel('L, cm');ylabel('T, K');grid on;

figure(4);plot(L*100, u);hold on;plot(L*100, u, 'ro');
title('Oxidizer Injection Analysis: Velocity');
xlabel('L, cm');ylabel('u, m/s');grid on;

figure(5);plot(L*100, M);hold on;plot(L*100, M, 'ro');
title('Oxidizer Injection Analysis: Mach Number');
xlabel('L, cm');ylabel('M');grid on;