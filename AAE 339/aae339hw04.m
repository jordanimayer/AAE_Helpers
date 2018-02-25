% AAE 339 - HW 04
% Jordan Mayer

% Problem 2

negHs = zeros(1,1000);
k=0;
for t=linspace(1,500,1000)
    k=k+1;
    negH = 0.00242*t^3 - 1.3548*t^2 + 6.32667*t^1.5 - 58.72*t^1.25;
    negH = negH + 202.573*t + 3845.92*t^(-0.5) - 4022.63*t^(-1);
    negH = negH + 1538.2*t^(-2) - 1445.56;
    negH = negH + t^(-2)*(-3.2524*10^(-14)*t^(2.5)) - 6.5048*10^(-14)*log(t);
    negHs(k) = negH;
    if (abs(negH - 2044) < 1)
        disp(negH);
        disp(t);
        break;
    end
end
temp=linspace(1,1000,1000)*100;
plot(temp,negHs);