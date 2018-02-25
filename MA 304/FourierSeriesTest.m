x = linspace(0, 2, 41)';
y = mod(2*x, 1);
[a, b, yfit] = Fseries(x, y, 10);
xfine = linspace(0, 2)';
yfine = Fseriesval(a, b, xfine);
plot(x, y, 'x', x, yfit, 'o', xfine, yfine);