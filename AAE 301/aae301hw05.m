% AAE 301 - HW 05
% Jordan Mayer
% Problem 8

b = [17, 176, 1236, 3356, 5239];
a = [1, 19, 188, 836, 1771, 1105];

[r, p, k] = residue(b,a)
t=linspace(0,7,1000);
g=zeros(1,1000);
for k=1:1000
    g(k)=2*exp(-6*t(k))*(4*cos(7*t(k))-5*sin(7*t(k)));
    g(k)=g(k)+2*exp(-3*t(k))*(2*cos(2*t(k))-3*sin(2*t(k)));
    g(k)=g(k)+5*exp(-1*t(k));
end
impulse(b,a,'^c');grid;hold on;
plot(t,g,'k');title('Calculated g(t) (pg. 109 #8)');
xlim([0 7]);xlabel('t');ylabel('g(t)');legend('using impulse','using residue');
hold off;