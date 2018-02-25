% AAE 333: HW 2 - Problem 3
% Plot C_D vs Re for a given sphere with unknown V_infty

rho_infty = 1.226;      %freestream density, kg/m^3
D = 3.70*10^(-4);       %drag force, N
S = 3.14*10^(-6);       %reference area, m^2
mu_infty = 1.780*10^(-5);  %freestream viscosity, kg/(m*s)
d_1 = 0.002;            %diameter of first sphere, m
d_2 = 0.003;            %diameter of second sphere, m

C_D_1 = transpose(linspace(0.06, 300, 1000));   %Possible drag coefficients for first sphere, dimensionless
C_D_2 = transpose(linspace(0.06, 300, 1000));   %Possible drag coefficients for second sphere, dimensionless
Re_1 = ones(length(C_D_1), 1) .* ((rho_infty * sqrt(2*D / (rho_infty*S)) * d_1) / mu_infty);   %Possible Reynolds numbers for first sphere, dimensionless
Re_1 = Re_1 ./ (C_D_1 .^ 0.5);
Re_2 = ones(length(C_D_2), 1) .* ((rho_infty * sqrt(2*D / (rho_infty*S)) * d_2) / mu_infty);   %Possible Reynolds numbers for second sphere, dimensionless
Re_2 = Re_2 ./ C_D_2 .^ 0.5;

loglog(Re_1, C_D_1);
hold on;
loglog(Re_2, C_D_2, '-ro');
axis([10^(-1), 10^7, 0.06, 400])
xlabel('Reynolds Number, dimensionless');
ylabel('Coefficient of Drag, dimensionless');
title('Drag on Falling Spheres');
legend('Diameter = 0.002 m', 'Diameter = 0.003 m');
grid on;