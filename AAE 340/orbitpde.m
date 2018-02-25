% evaluate state space system of differential equations
% described in 2b

function f= orbitpde (t,y)

    mu=3.986004418e14;
    y1 = y(1);
    y2 = y(2);
    y3 = y(3);
    y4 = y(4);

    y_dot1 = y2;
    y_dot2 = y1*y4^2 - mu*y1^(-2);
    y_dot3 = y4;
    y_dot4 = -2*y1^(-1)*y2*y4;
    
    f=[y_dot1 y_dot2 y_dot3 y_dot4]';
    
end