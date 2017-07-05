syms x1 x2 y1 y2 c

poly = expand((x1*y1 + x2*y2 + c)^2);
cy1 = coeffs(poly, [y1,y2])

