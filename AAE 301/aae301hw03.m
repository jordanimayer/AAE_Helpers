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
% figure(1);
% plot(f_i);grid;pbaspect([1 1 1]);
% 
% figure(2);
% plot(f_ii);grid;pbaspect([1 1 1]);
% 
% figure(3);
% plot(f_iii);grid;pbaspect([1 1 1]);
% 
% figure(4);
% plot(f_iv);grid;pbaspect([1 1 1]);
% 
% figure(5);
% plot(f_v);grid;pbaspect([1 1 1]);
% 
% figure(6);
% comet(real(f_iv), imag(f_iv));grid;pbaspect([1 1 1]);
% 
% figure(7);
% comet(real(fft([0;0;3i], 4096)),imag(fft([0;0;-3i],4096)));grid;pbaspect([1 1 1]);
% 
% 
% 
% % Page 63
% % Problem 1
% 
% % Part a
% t = (0:2^14-1)*2*pi/(2^14);
% f=t.*(t<=pi/2)+(3*pi/2-2*t).*(pi/2<t)-(3*pi/2-2*t).*(pi<t)+(-5*pi/2+2*t).*(pi<t)-(-5*pi/2+2*t).*(3*pi/2<t)+(2*pi-t).*(3*pi/2<t);
% a=ifft(f);
% bar((0:10),abs(a(1:11)).^2);grid; %for checking only
% %largest frequencies: 0, 1, 2, 3, 5, 6, 10
% fprintf('a_0 = ');disp(a(1));
% fprintf('a_1 = ');disp(a(2));
% fprintf('a_2 = ');disp(a(3));
% fprintf('a_3 = ');disp(a(4));
% fprintf('a_5 = ');disp(a(6));
% fprintf('a_6 = ');disp(a(7));
% fprintf('a_10 = ');disp(a(11));
% 
% %part b
% 
% %part c
% p=zeros(1,2^14);
% p=p+real(a(1));
% for k=[1 2 3 5 6 10];p=p+2*real(a(k+1)*exp(-1i*k*t));end
% plot(t,f);grid;hold on;plot(t,p,'r'); %for checking only
% err = norm(ifft(f-p))^2
% 
% %part d
% b = zeros(1, 20); for k=1:20;b(k)=a(k+1);end
% b=[conj(b(20:-1:1)),a(1),b];
% hold off;bar((-20:20),abs(b).^2);grid;
% title('Power Spectrum for f (pg. 63 #1)');
% xlabel('k');
% ylabel('Power');
% 
% 
% 
% % Page 64
% % Problem 2
% 
% %part a
% t=(0:2^14-1)*2*pi/(2^14);
% f=sin(3.*t)+cos(sin(8.*t));
% a=ifft(f);
% bar((0:100),abs(a(1:101)).^2);grid; %for checking only
% %significant frequencies: 0, 3, 16
% fprintf('c_0 = ');disp(a(1));
% fprintf('c_3 = ');disp(a(4));
% fprintf('c_16 = ');disp(a(17));
% 
% %part b
% 
% %part c
% p=zeros(1,2^14);
% p=p+real(a(1));
% p=p+2*real(a(4)*exp(-1i*3*t));
% p=p+2*real(a(17)*exp(-1i*16*t));
% plot(t,f);grid;hold on;plot(t,p,'r'); %for checking only
% err = norm(ifft(f-p))^2
% 
% %part d
% b = zeros(1, 20); for k=1:20;b(k)=a(k+1);end
% b=[conj(b(20:-1:1)),a(1),b];
% hold off;bar((-20:20),abs(b).^2);grid;
% title('Power Spectrum for f (pg. 64 #2)');
% xlabel('k');
% ylabel('Power');
% 


% Page 64
% Problem 3

%part a
t=(0:2^14-1)*2*pi/(2^14);
f=sin(40*t).*exp(sin(8*t));
a=ifft(f);
bar((20:60),abs(a(21:61)).^2);grid; %for checking only
%significant frequencies: 24, 32, 40, 48, 56
fprintf('c_24 = ');disp(a(25));
fprintf('c_32 = ');disp(a(33));
fprintf('c_40 = ');disp(a(41));
fprintf('c_48 = ');disp(a(49));
fprintf('c_56 = ');disp(a(57));

% %part b
p=zeros(1,2^14);
p=p+2*real(a(25)*exp(-1i*24*t));
p=p+2*real(a(33)*exp(-1i*32*t));
p=p+2*real(a(41)*exp(-1i*40*t));
p=p+2*real(a(49)*exp(-1i*48*t));
p=p+2*real(a(57)*exp(-1i*56*t));
plot(t,f);grid;hold on;plot(t,p,'r');
err = norm(ifft(f-p))^2
 
%part c
b = zeros(1, 60); for k=1:60;b(k)=a(k+1);end
b=[conj(b(60:-1:1)),a(1),b];
hold off;bar((-60:60),abs(b).^2);grid;
title('Power Spectrum for f (pg. 64 #3)');
xlabel('k');
ylabel('Power');


% 
% %Page 64
% %Problem 4
% 
% t=(0:2^20-1)*1/(2^20);
% f=exp(t);
% g=t;
% a=ifft(f);
% b=ifft(g);
% 
% %part a
% int_f=norm(a)
% 
% %part b
% int_f2=norm(a)^2
% 
% %part c
% int_fg=a*b'