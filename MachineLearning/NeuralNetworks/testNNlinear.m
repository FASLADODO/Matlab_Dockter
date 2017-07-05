%test NN with linear data

clear all

%% make fake data 
nn = 200;
noiz = 0.7;

rangen = [0 5];

params1 = [4 15]; %m, b
params2 = [2 10]; %m, b


xm = linspace(rangen(1),rangen(2),nn)';

xm1 = xm + randn(nn,1)*noiz;
xm2 = xm + randn(nn,1)*noiz;

yx1 = xm1.*params1(1) + params1(2) + randn(nn,1)*noiz;
yx2 = xm2.*params2(1) + params2(2) + randn(nn,1)*noiz;


X1 = [xm1, yx1];
X2 = [xm2, yx2];

X = [X1;X2];
Y = [ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];
classes = unique(Y)


figure
gscatter(X(:,1),X(:,2),Y)
title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast');
set(hl,'FontSize',12);
title('data wtih true class')

%% see what the scaled data is
[Xhat,~,~] = MeanVarianceScale(X);
figure
gscatter(Xhat(:,1),Xhat(:,2),Y)
title('scaled data')

%%

% format output
[YBinary,cslist,mapping] = NNFormatOutput(Y);

%initilize the layers / weights
layers = 3;
nodes = [4,5,3];
[NN] = NNInitialize(X,YBinary,layers,nodes);

alpha = 0.99999; %learning rate (0 - 1)
max_epoch = 1000;
accuracy_target = 0.5;

%Train it 
[NN,MSE] = NNTrain(X,YBinary,NN,alpha,max_epoch,accuracy_target);


%% Check classifiying

%online NN classify
Yestb = NNOnline(NN, X);

%convert back to our class labels
Yest = NNUnformatOutput(Yestb,cslist,mapping);
    
figure
gscatter(X(:,1),X(:,2),Yest)
title('data wtih est class')













