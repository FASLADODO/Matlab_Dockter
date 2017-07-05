function K=gaussian_kern(xs,x,h)

% Gaussian kernel function
K=eusdist(diag(1./h)*x',diag(1./h)*xs');
K=exp(-K/2);