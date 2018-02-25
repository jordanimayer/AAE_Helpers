%%%%%
% AAE 364
% HW 06
% Jordan Mayer
%%%%%

%%% Problem 1. B-5-27 %%%

% initialize variables
J = 10; % assume constant J (only analyzing trends; no exact value needed)
K = [1 10 100];     % small, medium, large K
B = [1 10 100];     % small, medium, large B
t = linspace(0,2*pi,100);  % time vector
r = t;              % unit ramp input
y = zeros(100,9);  % response vectors (in columns)

% calculate responses
for i=1:3
    for j=1:3
        k = (i-1)*3+j;
        num = [0 0 K(i)];   % CLTF numerator
        den = [J B(j) K(i)];    % CLTF denominator
        y(:,k) = lsim(num, den, r, t);   % response
    end
end

% create plots to view trends
plot(t, r, '-b');   % input
hold on;
plot(t, y(:,1), 'og');  % response at K=1, B=1
plot(t, y(:,2), '^g');  % response at K=1, B=10
plot(t, y(:,3), 'sg');  % response at K=1, B=100
plot(t, y(:,4), 'or');  % response at K=10, B=1
plot(t, y(:,5), '^r');  % response at K=10, B=10
plot(t, y(:,6), 'sr');  % response at K=10, B=100
plot(t, y(:,7), 'ok');  % response at K=100, B=1
plot(t, y(:,8), '^k');  % response at K=100, B=10
plot(t, y(:,9), 'sk');  % response at K=100, B=100
hold off;
title('Unit-ramp Response at Varying (K, B)');
xlabel('t');ylabel('y(t)');grid on;
legend('r(t)','(1,1)','(1,10)','(1,100)','(10,1)','(10,10)','(10,100)',...
    '(100,1)','(100,10)','(100,100)','Location','Best');
% NOTE: when viewing plots:
    % K is grouped by color, B is grouped by marker
    % (this is to help in assessing trends)
    % blue -> K = 1
    % red -> K = 10
    % black -> K = 100
    % circle -> B = 1
    % triangle -> B = 10
    % square -> B = 100