function W = localweights(X, x_p, T)
    W = exp(- ( ( X - x_p).^2)/(2*T^2) );
end