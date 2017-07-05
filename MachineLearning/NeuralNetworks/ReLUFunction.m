function [f] = ReLUFunction(x)
%compute rectified linear unit activation function
% range: 0<f<inf

f = x;
f(x < 0) = 0;

end