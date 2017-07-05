% test planar division classifier for boosting

%points
d1 = [1, 3.5;
    1,7;
    2.5,4;
    3.6,3.5;
    4,1.5;
    5,4];

d2 = [2.5,6;
    3.5,7.1;
    5,6;
    7,1;
    7.5,3];

labels = [ones(length(d1),1)*(-1); ones(length(d2),1)*1 ];
Data = [d1; d2];

%plot true class
scale = 1000;
figure
gscatter(Data(:,1), Data(:,2), labels)
axis([0 8 0 8])
title('class true')

%% rods adaboost

%Rods custom training method
RodTree = AdaBoostTrain(Data,labels,5,true);

%classify using the Tree
[Est,Value] = AdaBoostClassify(Data,RodTree);

figure
gscatter(Data(:,1), Data(:,2), Est)
title('class estimate')

labels
Est
Value


%% now rand data

nn = 100;

mu1 = [ 1,2];
mu2 = [ 4,5];
sigma1 = [1, 0.1; 0.1 ,2.1];
sigma2 = [1, 0.1; 0.1 ,1.5];

data1 = mvnrnd(mu1,sigma1,nn);
data2 = mvnrnd(mu2,sigma2,nn);

Data = [data1;data2]; %2D data

labels = [ones(nn,1)*(-1); ones(nn,1)*1 ];

%plot true class
scale = 1000;
figure
gscatter(Data(:,1), Data(:,2), labels)
title('true class')

%% rods adaboost

rounds = 5;

%Rods custom training method
Tree = AdaBoostTrain(Data,labels,rounds,false);

%classify using the Tree
[Est,Value] = AdaBoostClassify(Data,Tree);

figure
gscatter(Data(:,1), Data(:,2), Est)
title('Est class')

corr = Est == labels;
acc = sum(corr)/length(corr)


%% Try it on circular data

theta = [0:0.01:pi/2]';
nn = length(theta);

r = 3.4;
x1 = [r*sin(theta) + randn(nn,1)*0.5, r*cos(theta) + randn(nn,1)*0.5];
x2 = [randn(nn,1),randn(nn,1)];

X = [x1;x2];
Labels = [ones(nn,1)*-1;ones(nn,1)*1];

figure
gscatter(X(:,1),X(:,2),Labels)

[Tree] = AdaBoostTrain(X,Labels,50,true);

[ClassEstimate,Value] = AdaBoostClassify(X,Tree);

figure
gscatter(X(:,1),X(:,2),ClassEstimate)
title('class est data')


