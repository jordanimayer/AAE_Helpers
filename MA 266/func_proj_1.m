function [ xp ] = func_proj_2( t, x )

xp = zeros(2,1);    %output is column vector
xp(1) = x(2);
xp(2) = -4*x(2)-5*x(1)+10*cos(0*t);

end

