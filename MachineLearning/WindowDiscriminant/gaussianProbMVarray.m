function P = gaussianProbMVarray(X,Sigma,Mu)

    sigin = inv(Sigma);
    ss = length(Sigma);
    [nn,~] = size(X);
    P = (1/ ( (2*pi)^(ss/2) .*sqrt(norm(Sigma)) ) ).* exp(-(1/2).* sum((X - repmat(Mu,nn,1))' .* (sigin*(X - repmat(Mu,nn,1))'), 1));
    
    P = P';
end