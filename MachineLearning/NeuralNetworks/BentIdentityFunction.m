function [f] = BentIdentityFunction(x)
%compute bent identity activation (-inf<x<inf)

f = ( ( sqrt(x.^2 +1) - 1) ./ 2) + x;

end