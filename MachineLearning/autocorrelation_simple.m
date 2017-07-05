function [R,good] = autocorrelation_simple(X,k)
%https://en.wikipedia.org/wiki/Autocorrelation
%see estimation section for discrete processes
%k must be less than length of n
    
    %key data info
    [NN,SS] = size(X);
    
    if(k >= NN)
       R = 0;
       good = false;
       return;
    end
    
    sigma = mean(std(X));
    mu = mean(X);
    
    %get X_t and X_t+k
    X1 = X(1:NN-k,:);
    X2 = X(k+1:NN,:);
    
    %zero mean shift
    shifted = bsxfun(@minus,X1,mu);
    stepshifted = bsxfun(@minus,X2,mu);
    
    %sum it only once to keep series
    %sumR = sum(sum(shifted.*stepshifted,2));
    sumR = sum(shifted.*stepshifted,2);
    
    %scale it
    R = (1 / ( (NN - k)*(sigma*2) ) ) .* sumR;
    
    good = true;
end