function P_bar = BayesianRegressionOnline(X,Y,Model)
%compute P_Ybar

[NN,~] = size(X);

for ii = 1:NN
    %row wise sample
    xbar = X(ii,:);
    ybar = Y(ii,:);
    
    %compute estimate
    ybar = Y(ii,:);

    %first we compute the mean estimate
    mu_bar = xbar*Model.mu_B;

    %then we compute the covariance of the sample
    sigma_bar = (Model.sigma.^2) + xbar*Model.sigma_B*xbar';

    %posterior probability of ybar given xbar
    p_ybar = NormDistND(ybar,mu_bar,sigma_bar);
    
    P_bar(ii,:) = p_ybar;
end

end