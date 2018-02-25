% % Page 58
% % Problem 3
% % Parts i-v, ix, x
% 
% f_i = fft([0;5], 4096);
% f_ii = real(f_i);
% f_iii = imag(f_i);
% f_iv = fft([0;0;0;-2i], 4096);
% f_v = conj(f_iv);
% 
% % figure(1);
% % plot(f_i);grid;pbaspect([1 1 1]);
% % 
% % figure(2);
% % plot(f_ii);grid;pbaspect([1 1 1]);
% % 
% % figure(3);
% % plot(f_iii);grid;pbaspect([1 1 1]);
% % 
% % figure(4);
% % plot(f_iv);grid;pbaspect([1 1 1]);
% % 
% % figure(5);
% % plot(f_v);grid;pbaspect([1 1 1]);
% 
% % figure(6);
% % comet(real(f_iv), imag(f_iv));grid;pbaspect([1 1 1]);
% 
% figure(7);
% comet(real(fft([0;0;3i], 4096)),imag(fft([0;0;-3i],4096)));grid;pbaspect([1 1 1]);


% Page 63
% Problem 1

% Part a
t = (0:2^14-1)*2*pi/(2^14);
f=t.*(t<=pi/2)+(3*pi/2-2*t).*(pi/2<t)-(3*pi/2-2*t).*(pi<t)+(-5*pi/2+2*t).*(pi<t)-(-5*pi/2+2*t).*(3*pi/2<t)+(2*pi-t).*(3*pi/2<t);
a=ifft(f);
bar((0:10),abs(a(1:11)).^2);grid;
%largest frequencies: 0, 1, 2, 3, 5, 6, 10
fprintf('a_0 = ');disp(a(1));
fprintf('2a_1 = ');disp(2*a(2));
fprintf('2a_2 = ');disp(2*a(3));
fprintf('2a_3 = ');disp(2*a(4));
fprintf('2a_5 = ');disp(2*a(6));
fprintf('2a_6 = ');disp(2*a(7));
fprintf('2a_10 = ');disp(2*a(11));

%part b

%part c
p=zeros(1,2^14);
p=p+real(a(1));
for k=1:10;p=p+2*real(a(k+1)*exp(-1i*k*t));end
plot(t,f);grid;hold on;plot(t,p,'r');
err = norm(ifft(f-p))^2

%part d
b = zeros(1, 20); for k=1:20;b(k)=a(k+1);end
b=[conj(b(20:-1:1)),a(1),b];
hold off;bar((-20:20),abs(b).^2);grid;
title('Power Spectrum for f (pg. 63 #1)');
xlabel('k');
ylabel('Power');