function g = sigmoid(X,theta)
    z = X*theta;
    g = 1/(1+exp(-z));
end