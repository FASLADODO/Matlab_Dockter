function f = gammaPDF(x,a,b)
%compute the gamma probability
%a and b are shape parameters

f = ( 1 / ( (b^a) * gamma(a) ) ) .* x.^(a-1) .* exp(-x./b);

end