%TestModelSelection

clear all

noiz = [0,0.2,1,1];
modelnoiz = 10;

nn = 100;
range = [1,35];

x0 = linspace(range(1),range(2),nn)';

params = [23;-4;-0.5;-0.1]

XAllog = [ones(nn,1), x0, x0.^2, x0.^3, x0.^4,  x0.^5,  x0.^6,  x0.^7];

Xreg = XAllog(:,[1:4]);

[NN,SS] = size(Xreg)

X = [];
Y = [];
for ii = 1:10
    rng('shuffle')
    Xt = Xreg + randn(NN,SS).*repmat(noiz,NN,1);
    Yt = Xt*params + randn(NN,1).*modelnoiz;
    X = [X; Xt];
    Y = [Y; Yt];
end


scatter(X(:,2),Y,'r.')

%%

[idxtrain, idxval, idxtest]  = dividerand(size(X,1),0.7,0.15,0.15);

%return sub sampled vector
trainX = XAll(idxtrain,:);
valX= XAll(idxval,:);
testX = XAll(idxtest,:);
trainY = Y(idxtrain,:);
valY= Y(idxval,:);
testY = Y(idxtest,:);

figure
scatter(trainX(:,2),trainY,'r.')
title('train')

figure
scatter(testX(:,2),testY,'b.')
title('validate')

%% functionized
rng('shuffle')

XAll = [X, X(:,2).^4,  X(:,2).^5,  X(:,2).^6,  X(:,2).^7];

[bestModel, minRMS, bestParams,bestCoords,~] = ModelSelection(XAll,Y)

Xro = XAllog(:,[bestModel]);
Yro = Xro * bestParams;

figure
plot(Xro(:,2),Yro,'r-+')
hold on
scatter(X(:,2),Y,'b.')
hold off
legend('best model','train')


