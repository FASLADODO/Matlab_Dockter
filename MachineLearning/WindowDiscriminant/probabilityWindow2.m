function [P] = probabilityWindow2(Data,CPoint)
    sigma = cov(Data);
    mu = mean(Data,1);
    
    P = gaussianScale(CPoint,sigma,mu);
end