%%%%%
% AAE 364
% HW 08
% Jordan Mayer
% 
% Perform various calculations and create plots for problems in HW 09.
% Output values (from disp) are written down in the hand-written portion
% of the assignment.
%%%%%

close all;  % close all open figures

%%% Problem 1: B-6-17 %%%
% Part a
% Use magnitude condition to find lead constant
sd = -2 + 1i*2*sqrt(3); % desired pole
Kl = abs(sd+4)*abs(sd)*abs(0.5*sd+1)/(5*abs(sd+2));	% lead constant
fprintf('\nProblem 1: B-6-17, Part a\nlead constant\n');
disp(Kl);

% Part b
% Plot unit step response
G = tf([0 0 5], [0.5 1 0]); % uncompensated transfer function
Gc = tf([1.6 3.2], [1 4]);  % compensator transfer function
sys_uc = feedback(G, 1);    % uncompensated control system
sys_c = feedback(G*Gc, 1);  % compensated control system
figure(1); step(sys_uc, '-k'); hold on; step(sys_c, '-r');
legend('uncompensated', 'compensated'); grid on;
title('Unit Step Response for Problem 1: B-6-17, Part b');


%%% Problem 1: B-6-19 %%%
% Use magnitude condition to find lag constant
sd = -2 + 1i*2*sqrt(3); % desired pole
KL = abs(sd+0.01)*abs(sd)*abs(sd+4)/(16*abs(sd+0.05));  % lag constant
fprintf('\nProblem 1: B-6-19, Part a\nlag constant\n');
disp(KL);


%%% Problem 1: B-6-21 %%%
% use magnitude condition to find lead constant
sd = -2 + 1i*2*sqrt(3); % desired pole
Kl = abs(sd+20.02)*abs(sd)*abs(sd+5)/10;
fprintf('\nProblem 1: B-6-21\nlead constant\n');
disp(Kl);
% use magnitude condition to find lag constant
KL = abs(sd)*abs(sd+0.01)*abs(sd+5)*abs(sd+20.02)/(336.36*abs(sd+0.1488));
fprintf('\nProblem 1: B-6-21\nlag constant\n');
disp(KL);


%%% Problem 1: B-6-22 %%%
% use magnitude condition to find lead constant
sd = -1.333+3.478*1i;   % desired pole
Kl = abs(sd+1.461)*abs(sd)*abs(sd+1)*abs(sd+2)/(abs(sd+1.333)*abs(2*sd+1));
fprintf('\nProblem 1: B-6-22\nlead constant\n');
disp(Kl);


%%% Problem 1: B-6-26 %%%
% part a
% plot root locus in MATLAB to compare to hand-drawn sketch
L = tf([0 0 0 2], [1 2 2 0]);   % L(s)
figure(2); rlocus(L);           % plot root locus
title('Root Locus for Problem 1: B-6-26, Part a');
axis([-4 4 -4 4]); pbaspect([1 1 1]);
    % square plot for easy comparison to hand-drawn sketch
% part b
% plot vector of angle 60 degrees to find intersection with root locus
hold on;
k = 0:1000;
plot(k*exp(deg2rad(180-60)*1i));
title('Problem 1: B-6-26, Part b');
% calculate L(sd) to find necessary a
sd = -0.5+1i*0.8662;
L = 2/(sd*(sd^2+2*sd+2));
fprintf('\nProblem 1: B-6-26\nL(sd)\n');
disp(L);


%%% Problem 2 %%%
% find angle of uncompensated transfer function for desired pole
sd = -2 + 2.728*1i; % desired pole
G = (1.1057*sd+0.1900)/(sd^3+0.7385*sd^2+0.9009*sd);
fprintf('\nProblem 2\nuncompensated TF angle\n');
disp(rad2deg(angle(G)));
% use magnitude condition to find PD constant
Kpd = abs(sd^3+0.7385*sd^2+0.8008*sd)/(abs(sd+3.097)*abs(1.1057*sd+0.19001));
fprintf('\nProblem 2\nPD constant\n');
disp(Kpd);
% use magnitude condition to find PI constant
Kpi = abs(sd)/abs(sd+0.01);
fprintf('\nProblem 2\nPI constant\n');
disp(Kpi);