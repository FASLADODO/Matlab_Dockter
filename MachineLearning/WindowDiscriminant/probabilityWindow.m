function [P] = probabilityWindow(Data,CPoint)
    sigma = cov(Data);
    mu = mean(Data,1);
    
    P = gaussianProbMV(CPoint,sigma,mu);
end