%% make some 2d data

nn = 1000;

mu1 = 2;
sg1 = 0.2;
x1 = [randn(nn,1)*sg1 + mu1, randn(nn,1)*sg1 + mu1 ];

mu2 = 2.1;
sg2 = 0.8;
x2 = [randn(nn,1)*sg2 + mu2, randn(nn,1)*sg2 + mu2 ];

%combine it
Data = [x1; x2];
Labels = [ones(nn,1)*1; ones(nn,1)*2];

figure
gscatter(Data(:,1),Data(:,2),Labels)

%% Get single class rbfs and plot

CUSE = 1;
[Difference] = RelativeRBFSingleClass(Data,Labels,CUSE);

figure
gscatter(Data(:,1),Data(:,2),Labels)
hold on
Surface3D(Data(:,1),Data(:,2),Difference);
hold off
title('rbf class 1')


%% loop through all differences for rbf class 1 and see which one maximizes information gain

Direction = [1;2];

[Thresh,IG,Stash] = selectThresholdIG(Difference,Labels,Direction);
Thresh
IG

figure
scatter(Stash(:,1),Stash(:,2))
xlabel('threshold')
ylabel('IG')
title('threshold choice from information')

%%  estimate class and plot threshdold

[classest,sep] = RelativeRBFSingleClassOnline(Data,Labels,Data,Direction,Thresh);

zh = Thresh;
mu = 2;
sq = 1;


figure
gscatter(Data(:,1),Data(:,2),classest)
title('class estimates')

figure
gscatter(Data(:,1),Data(:,2),classest)
hold on
H = plotZplane([mu,mu],sq,zh);
hold on
Surface3D(Data(:,1),Data(:,2),sep);
hold off
title('IG Threshold')

corr = classest == Labels;

acc = mean(corr)




