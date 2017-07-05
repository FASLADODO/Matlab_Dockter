function P = gaussianScale(X,Sigma,Mu)

    sigin = inv(Sigma);
    [nn,~] = size(X);
    
    P = exp( -(1/2)*(X - Mu)*sigin*(X - Mu)' );
    
end