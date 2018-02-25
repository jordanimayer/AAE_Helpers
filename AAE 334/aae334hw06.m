% AAE 334
% HW 06
% Jordan Mayer
%
% Compute the coefficients of the Fourier series
% expansion of the circulation around a 3D
% symmetric airfoil section (alpha_ZL = 0)
% and use them to compute lift and induced drag.

% Wing characteristics
AR=8;   % Aspect Ratio
b = 16/3;   % wingspan
N = 100;
i = (linspace(1, N))';
theta = (i/(N+1))*pi;   % angle from left side (centered at mid-span)
y = -(b/2)*cos(theta);  % distance from mid-span
c = 1-(2*y/b).^2;   % chord along span

% Calculate CL, CDi, e for each wing
alpha_n = (pi/180)*linspace(-35, 35, N); % baseline angle-of-attack
CL_ut = zeros(a,1); %ut = untwisted
CL_tw = zeros(a,1); %tw = symmetrically twisted
CL_as = zeros(a,1); %as = asymmetrically twisted
CDi_ut = zeros(a,1);
CDi_tw = zeros(a,1);
CDi_as = zeros(a,1);
e_ut = zeros(a,1);
e_tw = zeros(a,1);
e_as = zeros(a,1);
for l=1:a
    % Create angle-of-attack matrices
    alpha_ut = ones(N,1) * alpha_n(l);
    alpha_tw = alpha_ut - (10*pi/180)*abs(2*y/b);
    alpha_as = alpha_ut - (10*pi/180)*(2*y/b);
    % Calculate NxN matrix
    X = zeros(N,N);
    for j=1:N
        for k=1:N
            X(j,k) = (2*b/(pi*c(j)) + k/sin(theta(j)))*sin(k*theta(j));
        end
    end
    % Calculate coefficients
    A_ut = inv(X)*alpha_ut;
    A_tw = inv(X)*alpha_tw;
    A_as = inv(X)*alpha_as;

    % Calculate CD_i, C_L, e
    sum_ut = 0;
    sum_tw = 0;
    sum_as = 0;
    for n=1:N
        sum_ut = sum_ut + n*A_ut(n)^2;
        sum_tw = sum_tw + n*A_tw(n)^2;
        sum_as = sum_as + n*A_as(n)^2;
    end
    CL_ut(l) = A_ut(1)*pi*AR;
    CL_tw(l) = A_tw(1)*pi*AR;
    CL_as(l) = A_as(1)*pi*AR;
    e_ut(l) = A_ut(1)^2/sum_ut;
    e_tw(l) = A_tw(1)^2/sum_tw;
    e_as(l) = A_as(1)^2/sum_as;
    CDi_ut(l) = CL_ut(l)^2/(pi*e_ut(l)*AR);
    CDi_tw(l) = CL_tw(l)^2/(pi*e_tw(l)*AR);
    CDi_as(l) = CL_as(l)^2/(pi*e_as(l)*AR);
end

% Part 1: plot CDi vs CL
figure(1);
plot(CL_ut, CDi_ut, 'b');hold on;
plot(CL_tw, CDi_tw, 'r');
plot(CL_as, CDi_as, 'k');
title('Elliptic Wing Analysis');
xlabel('CL');ylabel('CDi');xlim([-2 2]);grid on;
legend('Untwisted', 'Twisted', 'Asymmetric');

% Part 2: plot e vs CL
figure(2);
plot(CL_ut, e_ut, 'b');hold on;
plot(CL_tw, e_tw, 'r');
plot(CL_as, e_as, 'k');
title('Elliptic Wing Analysis');
xlabel('CL');ylabel('Efficiency Factor');xlim([-2 2]);
grid on;legend('Untwisted', 'Twisted', 'Asymmetric');

% Calculate circulation at alpha = 5 deg (for Part 3)

% Create angle-of-attack matrices
alpha_ut = ones(N,1) * 5*pi/180;
alpha_tw = alpha_ut - (10*pi/180)*abs(2*y/b);
alpha_as = alpha_ut - (10*pi/180)*(2*y/b);
% Calculate NxN matrix
X = zeros(N,N);
for j=1:N
    for k=1:N
        X(j,k) = (2*b/(pi*c(j)) + k/sin(theta(j)))*sin(k*theta(j));
    end
end
% Calculate coefficients
A_ut = inv(X)*alpha_ut;
A_tw = inv(X)*alpha_tw;
A_as = inv(X)*alpha_as;
% Calculate circulations
gamma_ut = zeros(N,1); %untwisted
gamma_tw = zeros(N,1); %symmetrically twisted
gamma_as = zeros(N,1); %asymmetrically twisted
gamma_el = zeros(N,1); %elliptical
for i=1:N
    for n=1:N
        gamma_ut(i) = gamma_ut(i) + A_ut(n)*sin(n*theta(i));
        gamma_tw(i) = gamma_tw(i) + A_tw(n)*sin(n*theta(i));
        gamma_as(i) = gamma_as(i) + A_as(n)*sin(n*theta(i));
    end
    gamma_el(i) = A_ut(1)*sin(theta(i));
end

% Part 3: plot circulation vs y/b 
figure(3);
plot(y/b,gamma_ut,'b');hold on;
plot(y/b,gamma_tw,'r');
plot(y/b,gamma_as,'k');
plot(y/b,gamma_el,'g');
title('Elliptic Wing Analysis');
legend('Untwisted','Twisted','Asymmetric','Eliptical');
xlabel('y/b');ylabel('Circulation');
grid on;