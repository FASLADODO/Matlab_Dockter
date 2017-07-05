function P = gaussianProbMV(X,Sigma,Mu)

    sigin = inv(Sigma);
    ss = length(Sigma);
    [nn,~] = size(X);
    
    scale = ( (2*pi)^(ss/2) *sqrt(norm(Sigma)) );
    
    P = (1/ ( (2*pi)^(ss/2) *sqrt(norm(Sigma)) ) ) * exp( -(1/2)*(X - Mu)*sigin*(X - Mu)' );
    
    P = P*scale;
end