function [Model] = BayesianRegressionTrain(X,Y)
%https://www.cs.utah.edu/~fletcher/cs6957/lectures/BayesianLinearRegression.pdf
%https://en.wikipedia.org/wiki/Bayesian_linear_regression

[~,k ] =size(X);

%OLS Solution
Model.B_OLS = pinv(X)*Y;

%now lets get the covariance
%first the estimate of y
mu_y = X*Model.B_OLS;

%residuals
R = Y - mu_y;

%covariance on data residuals
Model.sigma = std(R);

% a and b parameters for gamma distribution
a0 = 0.01;
Lambda_0 = a0 * ones(k);

%update parameters
%covariance on params
Lambda_n = X'*X + Lambda_0;

%mean on parameters
Model.mu_B = inv(Lambda_n)*X'*Y;

% covariance on parameters
Model.sigma_B = (Model.sigma.^2)*inv(Lambda_n);

%posterior estimate of B
Model.p_B = NormDistND(Model.B_OLS,Model.mu_B,Model.sigma_B );

end