function P = gaussianScale(Model)

    Sigma = Model.sigma;
    Mu = Model.mean;
    
    P = mvnpdf(Mu,Mu,Sigma);
end