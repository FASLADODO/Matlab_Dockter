%gaussian process for classification

%now we treat the probability as a sigmoid for -1/+1 class membership
%http://www.jmlr.org/papers/volume9/nickisch08a/nickisch08a.pdf

%The idea is we use relief-rbf to subsample d<n points from the training
%data with the highest separability (W_rbf)
%Then we use only those points in our classification

%% make up some data

% Now lets test gaussian process regression
fsize = 14;

%create some crazy data
nn = 1000;
noise = 0.5;

x1 = randn(nn,2) + repmat([2,2],nn,1);
x2 = randn(nn,2) + repmat([4,4],nn,1);
x3 = randn(nn,2) + repmat([6,6],nn,1);

X = [x1;x2;x3];
Y = [ones(nn,1)*-1;ones(nn,1)*1;ones(nn,1)*-1]; %labels

figure
gscatter(X(:,1),X(:,2),Y);
title('raw data')

%just using class membership as Y value
Model = GPRTrain(X,Y);
[Yest,Sigma,S2] = GPROnline(X,Model);
% gprMdl = fitrgp(X,Y);
% Yest = predict(gprMdl,X);

figure
scatter(X(:,1),X(:,2),10,Yest);
title('est class')
colorbar
colormap cool
%% we will use sigmoid logit funcion

siglogit = @(t) 1 ./ (1+exp(-t));

t = -10:0.1:10;

%change the latent function scale
sigma = 1;
%sigma = 0 -> sig(t) = 1/2 ignortant
%sigma = inf -> sig(t)=0; t<0, sig(t)=0.5; t=0; sig(t)=1; t > 0

figure
scatter(t,siglogit(sigma*t))

%% likelihood given sigmoid

%create data
nn = 1000;
x1 = randn(nn,2);
x2 = randn(nn,2) + repmat([2,2],nn,1);

X = [x1;x2];
Y = [ones(nn,1)*-1;ones(nn,1)*1];

%dummy function on data
f_x = NormRowWise(X-repmat([1,1],2*nn,1));

%likelihood
P_yx = siglogit(Y.*f_x);

figure
gscatter(X(:,1),X(:,2),Y);

figure
scatter(X(:,1),X(:,2),10,P_yx)
colorbar
colormap cool

%% lets try it?

%input args
trainratio = 0.5;
testratio = 0.2;

sigman = std(Y); %observation noise
sigmaf = mean(std(X)); %process noise
%silvermans rule of thumb
lscale = 1.06 * sigman * length(Y)^(-1/5) ; %length scale
kparams = [lscale, sigmaf];

%data size
[NN,SS] = size(X);

%first subsample our data
%all indices
idall = randperm(NN);

%get our train and test set
trainNN = floor(NN*trainratio);
testNN = floor(NN*testratio);
idtrain = idall( 1:trainNN );
idtest = idall( (trainNN+1):(trainNN+testNN) );

%grab a subset of the data to train with
XTrain = X(idtrain,:); %used in model creation
YTrain = Y(idtrain,:); %used in model creation
%for initial testing
XTest = X(idtest,:); %used to test model
[~,idsort] = sort(XTest(:,1)); %make our X values be sorted somehow
XTest = XTest(idsort,:);
YTest = Y(idtest(idsort),:); %not actually used

% now we compute our covariances, kernels, and means

%now lets get our 3 covariances (K_N,K_*N,K_**)
K_train_train = GaussianKernel(XTrain,XTrain,kparams);
K_test_train = GaussianKernel(XTest,XTrain,kparams);
K_test_test = GaussianKernel(XTest,XTest,kparams);

%noise on the covariance
Kns = K_train_train + sigman*eye(trainNN);


%WE STORE THIS FOR ONLINE SHIT
K_inv = inv(Kns);

%we compute this only once since we use it twice
% LK = K_test_train * K_inv;
LK = K_test_train/Kns;

%gaussian approximation to the mean and covariance
m_approx = YTrain;
cov_approx = K_inv;

%now we compute the mean
mu_test = K_test_train*K_inv*m_approx;
%posterior variance
kss = K_test_test + sigman*eye(testNN);
sigma_test = kss - (K_test_train(K_inv - K_inv*cov_approx*K_inv)*K_test_train);



%upper and lower bounds
yu = mu_test + stdv;
y1 = mu_test - stdv;




