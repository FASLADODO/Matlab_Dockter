function [f] = LeakyReLUDerivative(x)
%derivative of leaky rectified linear unit activation function
% range: -inf<f<inf

f = x;
f(x < 0) = 0.01;
f(x >= 0) = 1;

end