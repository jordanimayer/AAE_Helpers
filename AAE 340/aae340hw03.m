% AAE 340 - HW 03
% Jordan Mayer

%%% Problem 3 %%%

% Initialize variables
wn = 1.0;                   % natural frequency, rad/s
eps = 0.675;                % damping ratio, dimensionless
m = 4;                      % mass, kg
f0 = -1;                    % forcing amplitude, N
lamda0 = 2.0;               % initial displacement, m
lamdad0 = 1.0;              % initial velocity, m/s
w = [0.6 1.2 1.8];          % forcing frequency, rad/s
t = linspace(0,6*pi,1000);  % time, s
lamda = zeros(3,1000);      % displacement, m

for i=1:3               % for each omega
    for k=1:1000            % for each t
        % calculate steady-state response
        AF = 1/sqrt((1-(w(i)/wn)^2)^2+(2*eps*(w(i)/wn))^2);
        phi = atan(-2*eps*(w(i)/wn)/(1-(w(i)^2/wn^2)));
        sgn = 1;
        if wn-w(i) < 0
            sgn = -1;
        end
        A = AF*f0/(m*wn^2)*sgn;
        lamdas = A*cos(w(i)*t(k)+phi);
        
        % calculate transient response
        C1 = lamda0 - A*cos(phi);
        C2 = (lamdad0 + eps*wn*C1 + A*w(i)*sin(w(i)))/(wn*sqrt(1-eps^2));
        lamdat = C1*cos(wn*sqrt(1-eps^2)*t(k));
        lamdat = lamdat + C2*sin(wn*sqrt(1-eps^2)*t(k));
        lamdat = exp(-eps*wn*t(k))*lamdat;
        
        % calculate full response
        lamda(i,k) = lamdas + lamdat;
    end
    % display w, phi, AF, and w/wn
    fprintf('\nw = %f\nphi = %f\nAF = %f\nratio = %f\n', ...
        w(i), phi, AF, w(i)/wn);
end

% Create plots
close all;
plot(t,lamda(1,:),'r-');hold on;
plot(t,lamda(2,:),'k-');plot(t,lamda(3,:),'b-');
title('Forced Mass-Spring-Damper Response');
xlabel('Time, s');ylabel('Displacement, m');
legend('w = 0.6 rad/s', 'w = 1.2 rad/s', 'w = 1.8 rad/s');
grid on;hold off;