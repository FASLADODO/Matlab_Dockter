function p = gaussianprob(y, u, sigma)
    if(nargin == 2)
       sigma =1; 
    end
    p = (1/sqrt(2*pi)).*exp(-(1/2).*(y-u).^2);
end