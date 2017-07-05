function [K] = GaussianKernel(X1,X2,kparams)
%compute the gaussian kernel
%kparams = [l,sigma]
%l is length scale
%sigmaf is process noise

dist = pdist2(X1,X2,'squaredeuclidean');
K = kparams(2) .* exp( -dist ./ (2*kparams(1)^2));

end