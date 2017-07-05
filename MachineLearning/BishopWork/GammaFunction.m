function T = GammaFunction(x)
%jk this already exists in matlab gamma(x)

u = 1:10000000;

T = sum( (u.^(x-1)) .* exp(-u) );

end