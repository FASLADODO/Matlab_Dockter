function gd = sigmoidDerivative(t)
%derivative of sigmoid activation function
%http://math.stackexchange.com/questions/78575/derivative-of-sigmoid-function-sigma-x-frac11e-x
    g = 1./(1+exp(-t));
    gd = g.*(1-g);
end