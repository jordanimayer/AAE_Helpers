%%%
% AAE 364 - HW 05
% Jordan Mayer
%%%

%%% Problem 3c %%%
% For a y(t) obtained in 3b, plot y(t) versus t

% initialize vectors for y and t
t = linspace(-5,5,1000);
y = zeros(1000);

% calculate y(t) at each value of t
for k=1:1000
    if (t(k) < 0)
        y(k) = 0;
    else
        y(k) = -1 + 4*exp(-t(k)) - 3*exp(-2*t(k));
    end
end

% plot y(t) versus t and format
plot(t,y,'-k');grid on;
title('Output for Unit Step Input (Problem 3c)');
xlabel('t');ylabel('y(t)');