function [L] = likelihoodfunction(Y,X,Theta)
%http://www.robots.ox.ac.uk/~fwood/teaching/W4315_Fall2011/Lectures/lecture_3/lecture_3.pdf

    [NN,SS] = size(X);
    
    %estimate
    Y_bar = X*Theta;
    
    %sigma squares
    sigma2 = (1/NN) * sum( (Y-Y_bar).^2);
    
    %residual
    RSS = sum( (Y-Y_bar).^2);
    
    %likelihood function
    L = (1 / ( (2*pi*sigma2)^(NN/2) )) * exp( (-1/(2*sigma2)) * RSS );
end