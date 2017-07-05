function [f] = ReLUDerivative(x)
%compute rectified linear unit derivative function

f = x;
f(x < 0) = 0;
f(x >= 0) = 1;

end