function Q = QuadraticForm(X,Sigma)
    Q = sum( X' .* (Sigma*X'), 1);
end