%Page 73
%Problem 3

t=(0:2^11-1)'*5/(2^11);
g=3+6*cos(200*t)-8*sin(900*t)+10*randn(size(t));
figure(1);plot(t,g);grid;title('g(t) (Pg. 73 #3)');xlabel('t');ylabel('g(t)');
f=3+6*cos(200*t)-8*sin(900*t);
figure(2);plot(t,f);grid;title('f(t) (Pg. 73 #3)');xlabel('t');ylabel('f(t)');
a=ifft(g);
figure(3);bar((0:1000),abs(a(1:1001)).^2);grid;
title('Power Spectrum for g(t) (Pg. 73 #3)');xlabel('k');ylabel('Power');
%very high frequencies: 0, 159, 716
%amplitudes:    k=0     -->     10.3134
%               k=159   -->     8.0215
%               k=716   -->     13.0345
