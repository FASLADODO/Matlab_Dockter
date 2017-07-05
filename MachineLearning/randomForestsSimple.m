% test random forests
% http://www.listendata.com/2014/11/random-forest-with-r.html
% https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
% https://www.mathworks.com/help/stats/treebagger.html

%create rand data
dist = 2;
nn = 100;
d1 = [randn(nn,1), randn(nn,1), randn(nn,1)];
d2 = [randn(nn,1) + dist, randn(nn,1) + dist, randn(nn,1) + dist];

Data = [d1; d2];
Labels = [ones(nn,1)*1; ones(nn,1)*2];


figure
gscatter3(Data(:,1),Data(:,2),Data(:,3),Labels)
title('true class')

cols = {'x1';'x2';'x3'};


%% implementation with functions

%train forest
[Forest, oobavg] = RandomForestTrain(Data,Labels,cols,15,2);
oobavg

%try classifying
LabelEstimate = RandomForestOnline(Data,Forest);

figure
gscatter3(Data(:,1),Data(:,2),Data(:,3),LabelEstimate)
title('Est class')

corr = LabelEstimate == Labels;
acc = mean(corr)

%% plot first tree

%grab first tree and display for kicks
treeTemp = Forest.Tree{1};
treeTemp.oobError
figure
DecisionTreePlot(treeTemp,Labels);

