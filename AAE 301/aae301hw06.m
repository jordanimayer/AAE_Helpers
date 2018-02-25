% Jordan Mayer
% AAE 301 - HW 06

% Pg. 117

% Problem 1
num = [0 0 12/5];
den = [1 1/5 1];

[r, p, k] = residue(num, den)
t=linspace(0,60,1000);
y=zeros(1,1000);
for k=1:1000
    y(k)=2.412*exp(-0.1*t(k))*sin(0.995*t(k));
end

figure(4);
plot(t,y,'-r');hold on;
impulse(num, den,'^c');grid on;title('y(t) (Pg. 117, #1)');
xlabel('t');ylabel('y(t)');legend('y(t)', 'impulse');


% Pg. 121

% Problem 2
num = [0 0 60 0 5 -160];
den = [1 3.1 10.3 24.2 16 0];

[r, p, k] = residue(num, den);

figure(1);
impulse(num, den);grid on;title('y(t) (Pg. 121, #2)');
xlabel('t');ylabel('y(t)');

% Problem 3
num = [0 21284 0 3128 -1192 64];
den = [178542 309804519 3140881044 317032 8 0];

[r, p, k] = residue(num, den);

figure(2);
impulse(num, den);grid on;title('y(t) (Pg. 121, #3)');
xlabel('t');ylabel('y(t)');

% Problem 4
num = [0 4 0 8 0 64];
den = [1 13 -8 54 8 0];

figure(3);
impulse(num, den);grid on;title('y(t) (Pg. 122, #4)');
xlabel('t');ylabel('y(t)');

% Problem 6
num = [0 0 0 4];
den = [2 6 4 0];
t = linspace(0, 50, 1000);
y = zeros(0, 1000);
for k=1:1000
    y(k) = 1 - 2*exp(-t(k)) + exp(-2*t(k));
end
figure(5);
plot(t,y, '-r');hold on;
impulse(num, den, '^c');grid on;title('y(t) (Pg. 122, #6)');
xlabel('t');ylabel('y(t)');legend('y(t)','impulse');