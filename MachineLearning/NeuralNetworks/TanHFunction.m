function f = TanHFunction(x)
%basic Tanh activation function
%range: -1<f<1
    f = ( 2./ (1+exp(-2.*x)) )  - 1;
end