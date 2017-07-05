% test simple decision tree
% https://alliance.seas.upenn.edu/~cis520/wiki/index.php?n=Lectures.DecisionTrees

%create rand data
dist = 3;
nn = 100;
d1 = [randn(nn,1), randn(nn,1), randn(nn,1)];
d2 = [randn(nn,1) + dist, randn(nn,1) + dist, randn(nn,1) + dist];

Data = [d1; d2];
Labels = [ones(nn,1)*1; ones(nn,1)*2];


figure
gscatter3(Data(:,1),Data(:,2),Data(:,3),Labels)
title('true class')

%% build a decision tree

collabels = {'x1';'x2';'x3'};

tic
t = DecisionTreeTrain(Data,Labels,collabels);
toc

labelest = DecisionTreeOnline(Data,t);

figure
gscatter3(Data(:,1),Data(:,2),Data(:,3),labelest)
title('est class')

%% display the tree

DecisionTreePlot(t,Labels);




