function P = gaussianScaleArray(X,Sigma,Mu)

    sigin = inv(Sigma);
    [nn,~] = size(X);
    P = exp(-(1/2).* sum((X - repmat(Mu,nn,1))' .* (sigin*(X - repmat(Mu,nn,1))'), 1));
    
    P = P';
end