[t,x]=ode45('func_proj_1', [0, 20], [0, 1]);
plot(t, x(:,1));
xlabel('t');
ylabel('u(t)')
title('Solutions for E = -0.51');