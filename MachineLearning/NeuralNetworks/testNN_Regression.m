%% test NN regression


%make some data
%params
nn = 100;
noisex = 0.2;
noisey = 0.9;
m = [5;5];
b = [4];

%sim
x1 = linspace(5,20,nn)';
x2 = linspace(1,10,nn)';
X = [x1,x2] + randn(nn,2).*noisex;
y1r = X*m + repmat(b,nn,1);
Y = y1r + randn(nn,1)*noisey;

figure
scatter3(X(:,1),X(:,2),Y)
title('raw data')


%% Lets try training it?

layers = 2;
nodes = [3,3];

NN = NNRegressionInitialize(X,Y,layers,nodes);

alpha = 0.1;
max_epoch = 1000;
mse_target = 0.1;

[NN,MSE] = NNRegressionTrain(X,Y,NN,alpha,max_epoch,mse_target);

%% return training regression

Y_Estimate = NNRegressionOnline(NN, X);

figure
scatter3(X(:,1),X(:,2),Y_Estimate,'r')
hold on
scatter3(X(:,1),X(:,2),Y,'b')
hold off
title('regressions')

%% linear fit

minx = min(X);
maxx = max(X);
x1lin = linspace(minx(1),maxx(1),nn)';
x2lin = linspace(minx(2),maxx(2),nn)';
X_Linear = [x1lin,x2lin];

Y_Linear = NNRegressionOnline(NN, X_Linear);

figure
scatter3(X_Linear(:,1),X_Linear(:,2),Y_Linear,'r')
hold on
scatter3(X(:,1),X(:,2),Y,'b')
hold off
title('regressions')
