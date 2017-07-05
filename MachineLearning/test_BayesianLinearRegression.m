% test bayesian linear regression
%https://en.wikipedia.org/wiki/Bayesian_linear_regression


nn =100;
noise = 1.5;

x1 = linspace(0,10,nn)' + randn(nn,1)*noise;
X = [x1, ones(nn,1)];

%parameters mx+b
B_true = [3; 1.5];

Y = X*B_true + randn(nn,1)*noise;

figure
scatter(X(:,1),Y(:,1))
title('raw data')

%% Weaponized version

[Model] = BayesianRegressionTrain(X,Y);

P_bar =  BayesianRegressionOnline(X,Y,Model);

% plot color map of probs
figure
scatter(X(:,1),Y(:,1),10,P_bar)
colorbar
colormap cool
title('posterior predictive')

%% old version with code
%First lets compute the true coefficients

B = pinv(X)*Y;

%now lets get the covariance
%first the estimate of y
mu_y = X*B;

%residuals
R = Y - mu_y;

%covariance on data
sigma = std(R)

%% Therefore we can compute basic gaussian regression as

%probability (x-u)=(y-f(x))
P_fx =  normpdf(Y,mu_y,sigma);

% plot color map of probs
figure
scatter(X(:,1),Y(:,1),10,P_fx)
colorbar
colormap cool
title('linear gaussian')

%% now bayesian regression we will give posterior for B

%we can compute our uncertainty in our parameter B
%by taking the mean (mu_n) and covariance (sigma_n) 
%B ~ N(mu_n,sigma_n)

% a and b parameters for gamma distribution
a0 = 0.01;
b0 = 0.01;
n = length(X);
k = length(B);
mu_0 = zeros(k,1); %prior mean
Lambda_0 = a0 * ones(k);

%update parameters
%covariance on params
Lambda_n = X'*X + Lambda_0;

%mean on parameters
mu_n = inv(Lambda_n)*X'*Y;

% covariance on parameters
sigma_n = (sigma.^2)*inv(Lambda_n);

%posterior estimate of B
p_B = NormDistND(B,mu_n,sigma_n)


%% Posterior Predictive Distribution

%now given a new sample xbar, we would like to predist the dependent
%value ybar given our parameters and uncertainty
%ybar ~ N(xbar*mu_n,sigma2)



for ii = 1:nn
    %row wise sample
    xbar = X(ii,:);
    ybar = Y(ii,:);

    %first we compute the mean estimate
    mu_bar = xbar*mu_n;

    %then we compute the covariance of the sample
    sigma_bar = (sigma.^2) + xbar*sigma_n*xbar';

    %posterior probability of ybar given xbar
    p_ybar = NormDistND(ybar,mu_bar,sigma_bar);
    
    P_bar(ii,:) = p_ybar;
end


% plot color map of probs
figure
scatter(X(:,1),Y(:,1),10,P_bar)
colorbar
colormap cool
title('posterior predictive')

