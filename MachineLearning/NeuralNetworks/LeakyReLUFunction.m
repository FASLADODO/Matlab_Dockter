function [f] = LeakyReLUFunction(x)
%compute leaky rectified linear unit activation function
% range: -inf<f<inf

f = x;
id = find(x<0);
f(id) = 0.01*x(id);

end