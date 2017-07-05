function [P] = probabilityWindowarray(Data)
    sigma = cov(Data);
    mu = mean(Data,1);
    
    P = gaussianProbMVarray(Data,sigma,mu);
end