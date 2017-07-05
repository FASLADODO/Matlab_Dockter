% 1D
nn = 50;

x1 = randn(nn,1);
x2 = randn(nn,1) + 2;

X = [x1;x2];
Labels = [ones(nn,1)*-1;ones(nn,1)*1];
W = ones(2*nn,1)/(2*nn);
Direction = [-1;1];

figure
gscatter(W,X,Labels)
title('original data')

[Thresh,bestscore,IG] = DecisionStumpBasic(X,W,Labels);
Thresh
IG
classest = DecisionStumpOnline(X,Thresh,Direction);
corr = classest == Labels;
acc = mean(corr)
[W,alpha] = AdaBoostUpdate(W,Labels,classest);


%Information gain decision stump
[Thresh2,IG2] = DecisionStumpIG(X,Labels)
classest2 = DecisionStumpOnline(X,Thresh,Direction);
corr2 = classest2 == Labels;
acc2 = mean(corr2)

figure
scale = 2000;
DT = getClassData(X,classest);
WT = getClassData(W,classest);
cc = 1;
scatter(ones(length(DT{cc}),1), DT{cc}, WT{cc}*scale,'ro')
hold on
cc = 2;
scatter(ones(length(DT{cc}),1), DT{cc}, WT{cc}*scale,'bo')
hold on
plot([0.5,1.5],[Thresh,Thresh],'k')
hold off
title('weights + data')


figure
gscatter(ones(2*nn,1),X,classest)
hold on
plot([0.5,1.5],[Thresh,Thresh],'k')
hold off
title('estimated classes')

%% 1D with different variances to look at information gain

nn = 50;

x1 = randn(nn,1)*0.1+2;
x2 = randn(nn,1)*0.3+3.5;

X = [x1;x2];
Labels = [ones(nn,1)*-1;ones(nn,1)*1];
W = ones(2*nn,1)/(2*nn);
Direction = [-1;1];

figure
gscatter(W,X,Labels)
title('original data')

%standard decision stump
[Thresh,bestscore,IG] = DecisionStumpBasic(X,W,Labels);
Thresh
IG
classest = DecisionStumpOnline(X,Thresh,Direction);
[W,alpha] = AdaBoostUpdate(W,Labels,classest);


%Information gain decision stump
[Thresh2,IG2] = DecisionStumpIG(X,Labels)

figure
scale = 500;
DT = getClassData(X,classest);
WT = getClassData(W,classest);
cc = 1;
scatter(ones(length(DT{cc}),1), DT{cc}, WT{cc}*scale,'ro')
hold on
cc = 2;
scatter(ones(length(DT{cc}),1), DT{cc}, WT{cc}*scale,'bo')
hold on
plot([0.5,1.5],[Thresh,Thresh],'k')
hold off
title('weights + data')


figure
gscatter(ones(2*nn,1),X,classest)
hold on
plot([0.5,1.5],[Thresh,Thresh],'k')
hold off
title('estimated classes')

%% 2D
nn = 100;

x1 = [randn(nn,1), randn(nn,1)];
x2 = [randn(nn,1) + 2, randn(nn,1) + 2];

X = [x1;x2];
Labels = [ones(nn,1)*-1;ones(nn,1)*1];
W = ones(2*nn,1)/(2*nn);
Dir = [-1,1];

figure
gscatter(X(:,1),X(:,2),Labels)
title('original data')

[Thresh1,bestscore,IG] = DecisionStumpBasic(X(:,1),W,Labels);
classest = DecisionStumpOnline(X(:,1),Thresh1,Dir);
[W,alpha] = AdaBoostUpdate(W,Labels,classest);

[Thresh2,bestscore,IG] = DecisionStumpBasic(X(:,2),W,Labels);
classest = DecisionStumpOnline(X(:,2),Thresh2,Dir);
[W,alpha] = AdaBoostUpdate(W,Labels,classest);

figure
scale = 3000;
DT = getClassData(X,classest);
WT = getClassData(W,classest);
cc = 1;
scatter(DT{cc}(:,1), DT{cc}(:,2), WT{cc}*scale,'ro')
hold on
cc = 2;
scatter(DT{cc}(:,1), DT{cc}(:,2), WT{cc}*scale,'bo')
hold on
plot([-5,5],[Thresh1,Thresh1],'k')
hold on
plot([Thresh2,Thresh2],[-5,5],'k')
hold off
title('weights + data')

figure
gscatter(X(:,1),X(:,2),classest)
title('class est data')

%% Now try adaboost
nn = 100;

x1 = [randn(nn,1), randn(nn,1)];
x2 = [randn(nn,1) + 2, randn(nn,1) + 2];

X = [x1;x2];
Labels = [ones(nn,1)*-1;ones(nn,1)*1];

[Tree] = AdaBoostTrain(X,Labels,20,true);

[ClassEstimate,Value] = AdaBoostClassify(X,Tree);

figure
gscatter(X(:,1),X(:,2),Labels)
title('original data')

figure
gscatter(X(:,1),X(:,2),ClassEstimate)
title('class est data')

%% try it on fisher

load fisheriris.mat

[X,Labels] = GetPNLabels(meas,species,[2,3]);

figure
gscatter3(X(:,1),X(:,2),X(:,4),Labels)
title('-1/1  data ')

[Tree] = AdaBoostTrain(X,Labels,20,true);

[ClassEstimate,Value] = AdaBoostClassify(X,Tree);

figure
gscatter3(X(:,1),X(:,2),X(:,4),ClassEstimate)
title('class est data')

corr = ClassEstimate == Labels;
acc = mean(corr)









