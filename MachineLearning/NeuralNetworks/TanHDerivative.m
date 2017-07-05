function g = TanHDerivative(x)
%derivative of tanh activation function
    f = ( 2./ (1+exp(-2.*x)) )  - 1;
    g = 1 - f.^2;
end