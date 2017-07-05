function [f] = BentIdentityDerivative(x)
%compute bent identity derivative (-inf<x<inf)

f = ( x ./ ( 2.*sqrt(x.^2 +1) ) ) + 1;

end