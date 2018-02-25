% AAE 340 - HW 02
% Jordan Mayer

%%% Problem 2d %%%

% Calculate u(t)
t = linspace(0, 2*pi, 1000);
u = zeros(1,1000);
for k=1:1000
    u(k) = 0.33*cos(2.487*t(k)) - 0.01659*sin(2.487*t(k));
    u(k) = u(k)*exp(-0.25*t(k));
    u(k) = u(k) + 1.570;
end

% Create and format plot
figure(1);
plot(t, u, '-k');
title('Displacement for Mass-Spring-Damper System (Problem 2d)');
xlabel('t');ylabel('u(t)');grid on;


%%% Problem 3d %%%

% Calcuate lamda1, lamda2, lamda
t = linspace(0, 2*pi, 1000);
lamda1 = zeros(1,1000);
lamda2 = zeros(1,1000);
lamda = zeros(1,1000);
for k=1:1000
    lamda1(k) = 4.6426*exp(-0.5367*t(k));
    lamda2(k) = -2.1426*exp(-1.862*t(k));
    lamda(k) = lamda1(k) + lamda2(k);
end

% Create and format plot
figure(2);
plot(t, lamda1, '-r');hold on;
plot(t, lamda2, '-b');plot(t, lamda, '-k');
title('Response of Mass-Spring-Damper System (Problem 3d)');
xlabel('t');ylabel('response');grid on;
legend('first component','second component','net response');
hold off;