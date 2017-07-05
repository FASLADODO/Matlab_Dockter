function g = sigmoidFunction(t)
%basic sigmoid activation function
%range: 0<f<1
    g = 1./(1+exp(-t));
end