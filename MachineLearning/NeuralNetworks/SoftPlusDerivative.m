function [f] = SoftPlusDerivative(x)
%compute softplus derivative function

f = 1 ./ (1 + exp(-x)) ;

end