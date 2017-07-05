function [Y,Sigma,S2] = GPROnline(X,Model)
%train a gaussian process regression model using data
%https://pdfs.semanticscholar.org/9b21/6ab14159c9193d6b10ddb370c636e9557c67.pdf
%https://www.quora.com/How-can-I-use-Gaussian-processes-to-perform-regression
%https://en.wikipedia.org/wiki/Gaussian_process
%X test data matrix, rows as samples
%Model from GPRTrain.m

%data size
[NN,SS] = size(X);

%compute our covariances, kernels, and means
K_test_train = GaussianKernel(X,Model.XTrain,Model.kparams);
K_test_test = GaussianKernel(X,X,Model.kparams);

%Now lets get our mean and covariance for the prediction
%we compute this only once since we use it twice
LK = K_test_train * Model.K_inv;

%now we compute the mean
Y = LK * Model.YTrain ;
%posterior variance
Sigma = (K_test_test + Model.sigman.^2*eye(NN) ) - LK* K_test_train';


%For plots only
%lastly compute the standard deviation 
S2 = diag(K_test_test) - sum(LK.^2, 2);

end