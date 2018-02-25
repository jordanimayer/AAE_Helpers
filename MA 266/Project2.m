[t,x]=ode45('func_proj_2', [0, 80], [0, 0]);
plot(t, x(:,1));
xlabel('t');
ylabel('Q(t)')
title('Solutions for w = 2');