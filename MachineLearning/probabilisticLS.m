function P = probabilisticLS(X,Y,theta,sigma)
    
    expx = Y - X*theta;
    
    P = (1/sqrt(2*pi*sigma)) * exp( -(expx.^2)/(2*sigma^2) );

end