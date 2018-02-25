% % AAE 301 - Final Exam
% %
% % Jordan Mayer
% 
% Problem 1

% Find cos(omega*t)
omega = [2, -1; -1, 2];
syms t;
ex1 = expm(i*omega*t)
ex2 = expm(-i*omega*t)

% Find cos(omega*t)^2
exZero = expm(zeros(2,2))
ex1 = expm(2*i*omega*t)
ex2 = expm(-2*i*omega*t)


% Problem 2

% Find fourier expansion of u
syms t;
syms k;
T = 2*pi;
w_0 = 1;
u = t-pi;
a_k = int(u*exp(i*k*w_0*t),t,0,T);
a_k = a_k/T
a_0 = int(u,t,0,T);
a_0 = a_0/T
alpha = int((t-pi)*cos(k*w_0*t),t,0,T);
alpha = 2*alpha/T
beta = int((t-pi)*sin(k*w_0*t),t,0,T);
beta = 2*beta/T

% Plot u, y, yss
t = linspace(0,16*pi,2^10);
u = zeros(1,2^10);
yss = zeros(1,2^10); 
ag = zeros(1,2^10);
for k=1:2^10
    u(k) = t(k)-pi;
    s = k*1i;
    G = s/(5*s^2 + s + 5);
    ag(k) = angle(G);
    yss = yss + -2/k * abs(G) * sin(k*t - pi/2);
end
saw = pi*sawtooth(t);
num = [0, 0, 1];
den = [5, 1, 5];
figure(1);
lsim(num, den, saw, t); hold on;
plot(t,yss, '-r');grid on;hold off;
legend('u(t)', 'yss(t)');


% Problem 3

% Find G(s)
M = [2, -1; -1, 1];
phi = [2, 1; 1, 1];
K = [1, -1; -1, 2];
syms s;
b = [7; 3];
G = s^2 * M + s * phi + K;
G = inv(G);
G = G*b

% Find q1ss, q2ss
G1i = (7+4i)/(-4+13i);
G2i = (-3-7i)/(-4+13i);


% Problem 4
m1 = 1;
m2 = 2;
c1 = 2;
c2 = 3;
k1 = 6;
k2 = 8;
phi = [c1+c2, -c2, 0; -c2, c2, 0; 0, 0, 0];
syms s;

for m3=linspace(0.01, 100, 50)
    k3 = 4*m3;
    M = diag([m1, m2, m3]);
    K = [k1+k2, -k2, 0; -k2, k2+k3, -k3; 0, -k3, k3];
    A = M*s^2 + phi*s + K;
    d = det(A);
    p = poles(1/d,s);
    figure(2);
    hold on;
    if m3 == 0.01
        for i=1:6
            plot(p(i), 'xb');
        end
    else
        for i=1:6
            plot(p(i), '.b');
        end
    end
end
title('Poles of G(s) for varying m3');
xlabel('Real Part');ylabel('Imaginary Part');
xlim([-4, 1]);grid on;hold off;

m3 = 1.5;
k3 = 4*m3;
M = diag([m1, m2, m3]);
K = [k1+k2, -k2, 0; -k2, k2+k3, -k3; 0, -k3, k3];
A = M*s^2 + phi*s + K;
d = det(A);
p = poles(1/d,s)
hold on;
for i=1:6
    plot(p(i), 'xr');
    disp(p(i));
end
title('Poles of G(s) for varying values of m3');
xlabel('Real Part');ylabel('Imaginary Part');
xlim([-4, 1]);grid on;hold off;

t = linspace(0, 100, 2^20);
u = sin(2*t);
num = [0, 0, 0, m3*c2, m3*k2, k3*c2, k3*k2];
den = [70, 455, 6650, 27510, 76160, 4760, 6720];
figure(4);
lsim(num, den, u, t);grid on;hold off;
xlim([0, 20]);