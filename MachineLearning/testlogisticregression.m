%%Iris Data

load fisheriris

sp = nominal(species);

%getting numeric labels
classes = unique(species)
labs = 1:length(classes);
labsc = strread(num2str(labs),'%s')
labels = double(nominal(species,labsc));

%getting independent data
for ii = 1:length(classes)
    labelsall{ii} = labels(labels == ii);
    measall{ii} = meas(labels == ii,:);
end

%plot each class
figure
scatter3(measall{1}(:,1),measall{1}(:,2),measall{1}(:,3),'r.')
hold on
scatter3(measall{2}(:,1),measall{2}(:,2),measall{2}(:,3),'b.')
hold on
scatter3(measall{3}(:,1),measall{3}(:,2),measall{3}(:,3),'g.')
hold off


%train on just two classes
trainL = [labelsall{1}; labelsall{2}];
trainD = [measall{1}; measall{2}];

%get parameters
[ModelB,dev,stats] = mnrfit(trainD,trainL)

%get running data set
runD = [measall{2}; measall{1}];
runD = [ones(length(runD),1), runD];

%plot classification
classin1 = runD*ModelB(:,1); %class  1 vs 3 (setosa vs versicolor)

figure
plot(classin1)

%% Rand Data 2

%TestModelSelection
NN = 100;
SS = 2;

MU1 = [2,6];
MU2 = [3,2];
sig1 = 0.9;
sig2 = 0.9;

X1 = randn(NN,SS)*sig1 + repmat(MU1,NN,1);
X2 = randn(NN,SS)*sig2 + repmat(MU2,NN,1);

Data = [X1; X2];
Labels = [ones(NN,1)*0; ones(NN,1)*1];
% Labels = [repmat([1,0,0],NN,1); repmat([0,1,0],NN,1); repmat([0,0,1],NN,1) ];

figure
scatter(X1(:,1),X1(:,2),'r.')
hold on
scatter(X2(:,1),X2(:,2),'b.')
hold off

%Now Try and get params
ModelR =  LogRegTrain(Data,Labels);

[P1,L1] = LogRegOnline(Data,ModelR );

figure
scatter(Data(:,1),Data(:,2),20,P1(:,1))
colorbar
colormap cool


%% Rand Data 3

%TestModelSelection
NN = 100;
SS = 2;

MU1 = [2,6];
MU2 = [3,2];
MU3 = [6,0];
sig1 = 0.9;
sig2 = 0.9;
sig3 = 0.9;

X1 = randn(NN,SS)*sig1 + repmat(MU1,NN,1);
X2 = randn(NN,SS)*sig2 + repmat(MU2,NN,1);
X3 = randn(NN,SS)*sig3 + repmat(MU3,NN,1);

Data = [X1; X2; X3];
Labels = [ones(NN,1)*0; ones(NN,1)*1; ones(NN,1)*2];
% Labels = [repmat([1,0,0],NN,1); repmat([0,1,0],NN,1); repmat([0,0,1],NN,1) ];

figure
scatter(X1(:,1),X1(:,2),'r.')
hold on
scatter(X2(:,1),X2(:,2),'b.')
hold on
scatter(X3(:,1),X3(:,2),'g.')
hold off

%Now Try and get params
ModelR =  LogRegTrain(Data,Labels);

[P1,L1] = LogRegOnline(X1,ModelR );
[P2,L2] = LogRegOnline(X2,ModelR );
[P3,L3] = LogRegOnline(X3,ModelR );

figure
scatter(X1(:,1),X1(:,2),20,P1(:,1))
hold on
scatter(X2(:,1),X2(:,2),20,P2(:,2))
hold on
scatter(X3(:,1),X3(:,2),20,P3(:,3))
hold off
colorbar
colormap cool

%% Dynamic Data

%TestModelSelection

noiz = [0.01,0.1,0.01];
modelnoiz = 1;

nn = 1000;
range = [0,30];

x0 = linspace(range(1),range(2),nn)';

params1 = [8;2.8;0.9]
params2 = [8;2.4;0.7]

Xreg = [ones(nn,1), x0, x0.^2];

[NN,SS] = size(Xreg);

X1 = Xreg + randn(NN,SS).*repmat(noiz,NN,1);
X2 = Xreg + randn(NN,SS).*repmat(noiz,NN,1);

Y1 = X1*params1 + randn(NN,1).*modelnoiz;
Y2 = X2*params2 + randn(NN,1).*modelnoiz;

figure
scatter(X1(:,2),Y1,'r.')
hold on
scatter(X2(:,2),Y2,'b.')
hold off


%% I can do it myself asshole


Data = [X1 ,Y1; X2, Y2];
Labels = [ones(nn,1);ones(nn,1)*2];

%Now Try and get params
ModelB =  LogRegTrain(Data,Labels);

%do it online
DataOn = [X1 ,Y1];

[PALL1,Lm] = LogRegOnline(DataOn,ModelB );
[LLR,combos] = cumSumLLR(PALL1);

figure
scatter(1:length(Lm),Lm,'b.')
title('probabilities')
axis([0 1000 0 1])
title('logreg class 1')

figure
scatter(1:length(LLR),LLR,'r.')
title('LLR class 1')


%do it online
DataOn = [X2 ,Y2];

[PALL2,Lm] = LogRegOnline(DataOn,ModelB );
[LLR,combos] = cumSumLLR(PALL2);

figure
scatter(1:length(Lm),Lm,'b.')
title('probabilities')
axis([0 1000 0 1])
title('logreg class 2')

figure
scatter(1:length(LLR),LLR,'r.')
title('LLR class 2')

figure
scatter(1:length(PALL1(:,1)),PALL1(:,1),'r.')
hold on
scatter(1:length(PALL2(:,1)),PALL2(:,1),'g.')
hold off
title('probs')

%% Try over multiple runs

runs = 20;
noiz = [0.01,0.1,0.01];
modelnoiz = 1;
paramnoiz = 0.02;

nn = 100;
range = [0,30];

x0 = linspace(range(1),range(2),nn)';

params{1} = [8;2.8;0.9];
params{2} = [8;2.4;0.7];

Xreg = [ones(nn,1), x0, x0.^2];
[NN,SS] = size(Xreg);

DataRuns = [];
LabelRuns = [];
for cc = [1,2]
    for tt = 1:runs
        ParamsOn = params{cc} + randn(length(params{cc}),1)*paramnoiz;
        
        Xt = Xreg + randn(NN,SS).*repmat(noiz,NN,1);
        Yt = Xt*ParamsOn + randn(NN,1).*modelnoiz;

        DataRuns = [DataRuns; Xt, Yt];
        LabelRuns = [LabelRuns; repmat(cc,NN,1)];
    end
end

figure
gscatter(DataRuns(:,2), DataRuns(:,4), LabelRuns);

%% Do log reg for runs version

%Train params
ModelR =  LogRegTrain(DataRuns,LabelRuns);

%do it online
gc = 2;
DataOn = DataRuns(LabelRuns == gc,:);

[PALL,Lm] = LogRegOnline(DataOn,ModelR );
[LLR,combos] = cumSumLLR(PALL);

figure
scatter(1:length(Lm),Lm,'b.')
title('probabilities')
axis([0 1000 0 1])
title('logreg ')

figure
scatter(1:length(LLR),LLR,'r.')
title('LLR')





