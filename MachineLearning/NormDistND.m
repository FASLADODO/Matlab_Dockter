function p = NormDistND(X,mu,sigma)
%probability density at single ND point 

d = length(X);
p = (1 / (sqrt(norm(sigma) * (2*pi).^d)) ) * exp( -0.5* (X-mu)'*inv(sigma)*(X-mu));

end