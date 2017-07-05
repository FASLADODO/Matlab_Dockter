function f = InvGammaPDF(x,a,b)
%the inverse gamma function

f = ( (b^a) / gamma(a) ) .* ( x.^(-a-1) ) .* exp(-b./x);
end