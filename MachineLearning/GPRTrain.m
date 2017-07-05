function Model = GPRTrain(X,Y,trainratio,testratio,sigman,kparams)
%train a gaussian process regression model using data
%https://pdfs.semanticscholar.org/9b21/6ab14159c9193d6b10ddb370c636e9557c67.pdf
%https://www.quora.com/How-can-I-use-Gaussian-processes-to-perform-regression
%https://www.cs.toronto.edu/~hinton/csc2515/notes/gp_slides_fall08.pdf
%https://en.wikipedia.org/wiki/Gaussian_process

%X data matrix, rows as samples
%Y output column vector
%trainratio: 0.5 percent of points to use in training
%test ratio: 0.2 percent of points to test on
%sigman: 0.5 observation noise
%kparams: kernel parameters = [lengthscale, sigmaf] sigmaf: %process noise
%lengthscale: typical distance between peaks
%A conditional of a gaussian distribution is also gaussian

%check args
if(nargin < 4)
    trainratio = 0.5;
    testratio = 0.2;
end
if(nargin < 6)
    sigman = std(Y); %observation noise
    sigmaf = mean(std(X)); %process noise
    %silvermans rule of thumb
    lscale = 1.06 * sigman * length(Y)^(-1/5) ; %length scale
    kparams = [lscale, sigmaf];
end

Model = [];

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


Kns = K_train_train + sigman*eye(trainNN);

if(false)
    %We're gonna optimize our length scale
    %marginal likelihood (how good are our current params)
    plog_prev = -0.5*YTrain'*(Kns\YTrain) - 0.5*log10(det(Kns)) - (trainNN/2)*log10(2*pi);
    %refine length scale
    lprev = kparams(1);
    kparams(1) = lprev*((rand(1)*0.4)+0.8); %perturb length a little bit
    for oo = 1:100
        %recompute our covariance and likelihood
        K_train_train = GaussianKernel(XTrain,XTrain,kparams);
        Kns = K_train_train + sigman*eye(trainNN);
        plog = -0.5*YTrain'*(Kns\YTrain) - 0.5*log10(det(Kns)) - (trainNN/2)*log10(2*pi);

        %gradient ascent (maximize log likelihood)
        differp = plog - plog_prev;
        differk = kparams(1) - lprev;
        partial = differp / differk; %change in prob / change in params
        lprev = kparams(1);
        kparams(1) = lprev + partial * 0.001;

        %update plog
        plog_prev = plog;
    end
end

%Now lets get our mean and covariance for the prediction
%We could use cholesky decomposition (we use this for inverses)
%The inverse of a positive definite, hermitian matrix is just
%The inverses of its cholesky decompositions
%https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
% L = chol(K_train_train + sigman.^2*eye(trainNN),'lower');
% K_inv = inv(L)'*inv(L);

%WE STORE THIS FOR ONLINE SHIT
K_inv = inv(Kns);

%we compute this only once since we use it twice
% LK = K_test_train * K_inv;
LK = K_test_train/Kns;

%now we compute the mean
mu_test = LK * YTrain ;
%posterior variance
sigma_test = (K_test_test + sigman*eye(testNN) ) - LK* K_test_train';


%draw a sample from our posterior at our test points
% L2 = chol(K_test_test + sigman*eye(testNN) - LK*LK', 'lower');
% f_post = mu_test + L2* randn(testNN,1);

%For plots only
%lastly compute the standard deviation 
s2 = diag(K_test_test) - sum(LK.^2, 2);
stdv = sqrt(s2);

%upper and lower bounds
yu = mu_test + stdv;
y1 = mu_test - stdv;

%for plotting
Model.xshade = [XTest; flipud(XTest)];
Model.yshade = [yu; flipud(y1)];


%store everything in model
Model.K_inv = K_inv;
Model.XTrain = XTrain;
Model.YTrain = YTrain;
Model.kparams = kparams;
Model.sigman = sigman;
%not neccesary
Model.XTest = XTest;
Model.YTest = YTest;
%for plots
Model.mu_test = mu_test;
Model.sigma_test = sigma_test;
% Model.f_post = f_post;



end