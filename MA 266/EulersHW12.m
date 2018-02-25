[t, y] = eul('fcn1', [0, 2], 1, 1/22);
a = 0.0:0.5:3.0;
for i=1:1:length(t)
    l = t(i);
    if any(a==l)
        b=y(i);
        fprintf('%.5f\t%.5f\n', l, b)
    end
end