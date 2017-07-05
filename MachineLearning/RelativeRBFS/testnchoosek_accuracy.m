mu1 = [1,1,1,1];
sigma1 = [1,0.1,0.1,0.1;0.1,1,0.1,0.1;0.1,0.1,1,0.1;0.1,0.1,0.1,1];

mu2 = [1,3,4,5];
sigma2 = [1,0.1,0.1,0.1;0.1,1,0.1,0.1;0.1,0.1,1,0.1;0.1,0.1,0.1,1];

nn = 5000;

d1 = mvnrnd(mu1,sigma1,nn);
d2 = mvnrnd(mu2,sigma2,nn);

Data=[d1;d2];
Labels=[ones(nn,1)*1;ones(nn,1)*2];

figure
gscatter(Data(:,1),Data(:,2),Labels)

figure
gscatter(Data(:,2),Data(:,3),Labels)

figure
gscatter(Data(:,3),Data(:,4),Labels)

figure
gscatter3(Data(:,1),Data(:,3),Data(:,4),Labels)

%% SOMETHING WERID IS HAPPENING WITH SCALING

% TEST ACCURACY FOR EACH NUMBER OF STATES FIND BEST COMBO FOR EACH


X = Data(:,[3,4]);

[Difference,ClassData,ProbData] = SimpleRelativeRBFTrain(X,Labels);


figure
gscatter(X(:,1),X(:,2),Labels);
hold on
cc = 1;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),ProbData{cc}(:,1));
hold on
cc = 2;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),ProbData{cc}(:,1));
hold off
title('probs')

figure
gscatter(X(:,1),X(:,2),Labels);
hold on
cc = 1;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),Difference{cc});
hold on
cc = 2;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),Difference{cc});
hold off
title('seperability')

%%

collabz = {'1';'2';'3';'4'};
testcolumns = [1,2,3,4];

[Variations,BestSep] = RelativeRBFSeperability_Accuracy(Data,Labels,testcolumns,collabz,1,'ploton');
BestSep.bestcollabels
columnrbf = BestSep.bestcolumns;

