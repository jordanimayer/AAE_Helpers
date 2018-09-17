close all;

m = [2 1 1]'; l = [1 0.5]';  % P4
q_dot_0 = [0 0 0]';
q_0 = [0 -90.01 90]';  % IC4
q_0(2) = deg2rad(q_0(2));
q_0(3) = deg2rad(q_0(3));

[t, x] = sim('aae564hw01syssoln');
y = x(:,1);
theta1 = x(:,2);
theta2 = x(:,3);
figure();
plot(t,y); grid on; xlabel('t'); ylabel('y');
figure();
plot(t,theta1); grid on; xlabel('t'); ylabel('theta1');
figure();
plot(t,theta2); grid on; xlabel('t'); ylabel('theta2');