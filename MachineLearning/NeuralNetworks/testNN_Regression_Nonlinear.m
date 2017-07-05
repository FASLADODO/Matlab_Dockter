%% test NN regression


%make some data
%params
nn = 5000;
noisex = 0.01;
noisey = 0.05;
a = 4;
b = 0.1;

%sim
x1 = linspace(1,20,nn)';
X = [x1] + randn(nn,1).*noisex;
y1r = a.*exp(-b.*X);
Y = y1r + randn(nn,1)*noisey;

figure
scatter(X(:,1),Y)
title('raw data')


%% Lets try training it?

layers = 2;
nodes = [3,3];

%NN = NNRegressionInitialize(X,Y,layers,nodes);

alpha = 0.1;
max_epoch = 1000;
mse_target = 0.03;

[NN,MSE] = NNRegressionTrain(X,Y,NN,alpha,max_epoch,mse_target);

%% return training regression

Y_Estimate = NNRegressionOnline(NN, X);

figure
scatter(X(:,1),Y_Estimate,'r')
hold on
scatter(X(:,1),Y,'b')
hold off
title('regressions')

%% linear fit

minx = min(X);
maxx = max(X);
x1lin = linspace(minx(1),maxx(1),nn)';
X_Linear = [x1lin];

Y_Linear = NNRegressionOnline(NN, X_Linear);

figure
scatter(X_Linear(:,1),Y_Linear,'r')
hold on
scatter(X(:,1),Y,'b')
hold off
title('regressions')
