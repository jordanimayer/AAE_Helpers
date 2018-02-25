% AAE 340 - HW 01
% Jordan Mayer

% Problem 2c

% Verify hand-drawn plots

w = 1.5;                        % natural frequency, given
t = linspace(0,2*pi,100);       % time
x_i = zeros(1,100);             % x(t) for part i
x_ii = zeros(1,100);            % x(t) for part ii
x_iii = zeros(1,100);           % x(t) for part iii
% Calculate x(t) for parts i, ii, iii
for k=1:100
    x_i(k) = (2/1.5)*sin(1.5*t(k));
    x_ii(k) = cos(1.5*t(k));
    x_iii(k) = cos(1.5*t(k)) + (2/1.5)*sin(1.5*t(k));
end
% Create and format plots
figure(1);plot(t,x_i,'-k');
title('x(t) for Problem 2bi');xlabel('t');ylabel('x(t)');grid on;
figure(2);plot(t,x_ii,'-r');
title('x(t) for Problem 2bii');xlabel('t');ylabel('x(t)');grid on;
figure(3);plot(t,x_iii,'-b');
title('x(t) for Problem 2biii');xlabel('t');ylabel('x(t)');grid on;

% Verify equivalence of all solutions

w = 1.5;                        % natural frequency
x0 = 1;                         % initial displacement
xdot0 = 2;                      % initial velocity
x_form1 = zeros(1,100);         % x(t) using form 1
x_form2 = zeros(1,100);         % x(t) using form 2
x_form3 = zeros(1,100);         % x(t) using form 3
% Calculate x(t) using each form
for k=1:100
    x_form1(k) = x0*cos(w*t(k)) + xdot0/w*sin(w*t(k));
    x_form2(k) = sqrt(x0^2 + xdot0^2/w^2)*cos(w*t(k)+atan(-xdot0/(w*x0)));
    x_form3(k) = (x0/2 + xdot0/(2*1i*w))*exp(1i*w*t(k));
    x_form3(k) = x_form3(k) + (x0/2 - xdot0/(2*1i*w))*exp(-1i*w*t(k));
end
% Create and format plots in same figure
figure(4);plot(t,x_form1,'^c');hold on;
plot(t,x_form2,'+r');plot(t,x_form3,'-k');
title('Three EOM Solutions for Mass-Spring system');
xlabel('t');ylabel('x(t)');grid on;ylim([-2, 2]);
legend('Form 1', 'Form 2', 'Form 3');
hold off;


% Problem 3

% Verify hand-drawn plots

t = linspace(0,2*pi,100);       % time
x_i = zeros(1,100);             % x(t) for part i
x_ii = zeros(1,100);            % x(t) for part ii
% Calculate x(t) for parts i and ii
for k=1:100
    x_i(k) = 2.134*cos(1.061*t(k)-1.083);
    x_ii(k) = 1.539*cos(1.710*t(k)-0.8634);
end
% Create and format plots
figure(5);plot(t,x_i,'-b');
title('x(t) for Problem 3ai');xlabel('t');ylabel('x(t)');grid on;
figure(6);plot(t,x_ii,'-c');
title('x(t) for Problem 3aii');xlabel('t');ylabel('x(t)');grid on;