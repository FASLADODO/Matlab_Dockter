function [f] = SoftPlusFunction(x)
%compute softplus activation function
%range: 0<f<inf

f = log(1+exp(x));

end