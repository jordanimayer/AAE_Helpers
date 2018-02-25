% AAE 301 - HW 08
%
% Jordan Mayer

% Problem 5

num = [1, -2.0001];
den = [1, -1, -2];
g = tf(num, den);
figure(1);
step(1*g);grid on;title('Step response for u(t) = 1 (Pg. 168, #5, part a)');

figure(2);
num = [0, 1];
den = [1, 1];
g = tf(num, den);hold on;
step(1*g, '^c');
num = [0, 1, 2.01];
den = [1, 3, 2];
g = tf(num, den);
step(1*g, '-k');grid on;title('Step response for u(t) = 1 (Pg. 168, #5, part b)');
legend('exact', 'approximation');


% Problem 6
a = [-3, -6, -10, 10; 1, -1, 0, 0; 0, 1, 1, -2; 0, 0, 1, -2];
b = [1, 0, 0, 0]';
c = [1, -4, -2, 10];
d = 0;
[num, den] = ss2tf(a, b, c, d);g = tf(num, den)
g = minreal(g)
[A, B, C, D] = minreal(a, b, c, d)
figure(3);
step(30*g, '^c');hold on;

num = [0, 0, 30, -120, 90];
den = [1, 3, 7, 5, 0];
[r, p, k] = residue(num, den)
t = linspace(0, 12, 1000);
y = 42*exp(-t).*cos(2*t) + 6*exp(-t).*sin(2*t) - 60*exp(-t) + 18;
plot(t, y, '-k');title('y(t) for u(t)=30 (Pg. 168, #7, part g)');
xlabel('t');ylabel('y(t)');legend('step','y(t)');grid on;