function P = gaussianEval(X,Model)

    Sigma = Model.sigma;
    Mu = Model.mean;
    
    P = mvnpdf(X,Mu,Sigma);
    
    P = P ./ Model.scale;
end