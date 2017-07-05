function J = costfunction(Y,X,THETA)
    J = (1/2) * sum( (Y - X*THETA).^2 );
end