% Test gaussian process regression
%https://pdfs.semanticscholar.org/9b21/6ab14159c9193d6b10ddb370c636e9557c67.pdf
%https://www.quora.com/How-can-I-use-Gaussian-processes-to-perform-regression
%https://en.wikipedia.org/wiki/Gaussian_process

% built in matlab version


%make up some noisey data
n = 500;
x = linspace(0,2*pi,n)';

y = x.*sin(x) + 0.2*randn(n,1);

figure
scatter(x,y)


%fit a model
gprMdl = fitrgp(x,y,'Basis','linear',...
      'FitMethod','exact','PredictMethod','exact');
  
% estimate data
ypred = resubPredict(gprMdl);

figure
plot(x,y,'b.');
hold on;
plot(x,ypred,'r','LineWidth',1.5);
xlabel('x');
ylabel('y');
legend('Data','GPR predictions');
hold off


%% now lets try it ourselves

% Now lets test gaussian process regression
fsize = 14;

%create some crazy data
nn = 1000;
noise = 0.5;
X = linspace(0,9,nn)' + randn(nn,1)*noise;

b = 2;

Y = X.*sin(X) + b  + randn(nn,1)*noise;

figure
scatter(X,Y);
title('raw data')

%% weaponized version

%length scale: typical distance between peaks

trainratio = 0.5; %how much of data will be in our model
testratio = 0.2; %how much of data should we test with
kparams = [1,1.5]; %lengthscale,sigmaf
sigman = 2.5; %observation noise
Model = GPRTrain(X,Y,trainratio,testratio);

figure
fill(Model.xshade, Model.yshade, [.9 .9 .9], 'linestyle', 'none')
hold on
scatter(Model.XTrain,Model.YTrain,'ro');
hold on
scatter(Model.XTest,Model.mu_test,'b+');
hold off
title('GP mean and error bars')
xlabel('x','FontSize',fsize)
ylabel('y','FontSize',fsize)
hleg = legend('\sigma Bounds','D_{train}','\mu_{test}')
hleg.FontSize = fsize;

figure
scatter(Model.XTest,Model.f_post,'b+');
hold on
scatter(Model.XTrain,Model.YTrain,'ro');
hold off
title('GP test points posterior')


% redraw some new data
nn2 = 200;
Xon = linspace(0,9,nn2)' + randn(nn2,1)*noise;

b = 2;
Yon = Xon.*sin(Xon) + b  + randn(nn2,1)*noise;

%do it online
[Yest,Sigmaest,S2est] = GPROnline(Xon,Model);

RMSE = sqrt(mean( (Yon-Yest).^2))

figure
scatter(Xon,Yest,'b+');
hold on
scatter(Xon,Yon,10,S2est);
hold off
colormap cool
colorbar
title('GP compute new data draw')

%% lets test parameter estimation

theta_0 = [1,1];
ynn = length(Model.YTrain);
Ky = GaussianKernel(Model.YTrain,Model.YTrain,theta_0);
% log marginal likelihood
logp = -0.5*Model.YTrain' * inv(Ky + sigman*eye(ynn) ) * Model.YTrain - 0.5*log10(norm(Ky)) - (ynn/2)*log10(2*pi);


%% next lets try on gaussian noise
%like here:
%http://katbailey.github.io/post/gaussian-processes-for-dummies/

p2 = 1;
lng = sqrt( (1/p2) / 2);
kparams = [lng,0.5]; %lengthscale,sigmaf
sigman = 0.001; %observation noise

XAll = [];
nn = 100;
figure
for ii = 1:10
    xi = linspace(-10,10,nn)';
    K_ss = GaussianKernel(xi,xi,kparams);
    Ld = chol(K_ss + sigman*eye(nn),'lower');
    prior = Ld*randn(nn,1);
    plot(xi,prior)
    hold on
    XAll(:,ii) = [ xi ];
end
hold on


%now lets choose some anchor points
xtr = [-6;-4;-2;0;2;4;6];
ytr = [-2;0;3;-3;0;4;0]/2;
scatter(xtr,ytr,50,'ko')
hold off
title('sample from kernel')

XTrain = xtr;
YTrain = ytr;
YTest = [];

figure
scatter(xtr,ytr,50,'ko')
hold on
for ii = 1:10
    XTest = XAll(:,ii);
    
    %now lets do it
    K_train_train = GaussianKernel(XTrain,XTrain,kparams);
    K_train_test = GaussianKernel(XTrain,XTest,kparams);
    K_test_test = GaussianKernel(XTest,XTest,kparams);

    %Now lets get our mean and covariance for the prediction

    %first we get cholesky decomposition (we use this for inverses)
    %The inverse of a positive definite, hermitian matrix is just
    %The inverses of its cholesky decomposition
    %https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
    L = chol(K_train_train + sigman*eye(length(XTrain)),'lower');

    %WE STORE THIS FOR ONLINE SHIT
    %K_inv = inv(L'*L);
    LK = K_train_test \ L;
    Ly = L \ YTrain;
    Ls = L \ K_train_test; %inv(K_train_train)*K_train_test

    %now we compute the mean
    mu_test = LK * Ly;
    sigma_test = (K_test_test + sigman*eye(nn) ) - LK*Ls;
    
    s2 = diag(K_test_test) - sum(LK.^2, 2);
    stdv = sqrt(s2);
    
    % Draw samples from the posterior at our test points.
    L2 = chol(K_test_test + 0.00001*eye(nn),'lower');
    f_post = mu_test + L2*randn(nn,1);

     plot(XTest,mu_test);
     hold on
     plot(XTest, f_post);
     hold on
end
hold off


%% OLD SCRIPT VERSION


%input args
trainratio = 0.5; %how much of data will be in our model
testratio = 0.2; %how much of data should we test with
kparams = [1,1.5]; %lengthscale,sigmaf
sigman = 1.5; %observation noise

[NN,SS] = size(X);

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


%This if we want to computationally determine parameters
if(false)
sigman = std(YTrain); %observation noise
sigmaf = mean(std(XTrain)); %process noise
lscale = sqrt(mean(range(XTrain))); %length scale
kparams = [lscale, sigmaf]
end

figure
scatter(X,Y,'ro');
hold on
scatter(XTrain,YTrain,'c*');
hold on
scatter(XTest,YTest,'k+');
hold off
% title('training data')
xlabel('x','FontSize',fsize)
ylabel('y','FontSize',fsize)
hleg = legend('All Data','D_{train}','D_{test}')
hleg.FontSize = fsize;

%% now we compute our covarainces, kernels, and means

%now lets get our 3 covariances (K_N,K_*N,K_**)
K_train_train = GaussianKernel(XTrain,XTrain,kparams);
K_test_train = GaussianKernel(XTest,XTrain,kparams);
K_test_test = GaussianKernel(XTest,XTest,kparams);

%Now lets get our mean and covariance for the prediction

%first we get cholesky decomposition (we use this for inverses)
%The inverse of a positive definite, hermitian matrix is just
%The inverses of its cholesky decomposition
%https://makarandtapaswi.wordpress.com/2011/07/08/cholesky-decomposition-for-matrix-inversion/
L = chol(K_train_train + sigman.^2*eye(trainNN));
K_inv = inv(L'*L);
%we compute this only once since we use it twice
LK = K_test_train * K_inv;

%now we compute the mean
mu_test = LK * YTrain ;
%posterior variance
sigma_test = (K_test_test + sigman.^2*eye(testNN) ) - LK* K_test_train';


%lastly compute the standard deviation for plotting
s2 = diag(K_test_test) - sum(LK.^2, 2);
stdv = sqrt(s2);

yu = mu_test + stdv;
y1 = mu_test - stdv;
xshade = [XTest; flipud(XTest)];
yshade = [yu; flipud(y1)];

figure
fill(xshade, yshade, [.9 .9 .9], 'linestyle', 'none')
hold on
scatter(XTrain,YTrain,'ro');
hold on
scatter(XTest,mu_test,'b+');
hold off
% title('GP estimate')
xlabel('x','FontSize',fsize)
ylabel('y','FontSize',fsize)
hleg = legend('\sigma Bounds','D_{train}','\mu_{test}')
hleg.FontSize = fsize;









