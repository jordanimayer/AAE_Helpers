x = linspace(0, 10, 50);
c = linspace(0, 10, 20);

for i=1:1:20
    plot(x, c(i)./x);
    hold on;
    plot(x, sqrt(x.^2 - c(i)));
    %hold on;
    axis([0, 10, 0, 10]);
    title('Potential Flow in First Quadrant');
    xlabel('x (m)');
    ylabel('y (m)');
end