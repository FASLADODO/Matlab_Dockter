mu1 = [1,1,1];
sigma1 = [1,0.1,0.1;0.1,1,0.1;0.1,0.1,1];

mu2 = [4.5,4.5,4.5];
sigma2 = [1,0.1,0.1;0.1,1,0.1;0.1,0.1,1];

nn = 2000;

d1 = mvnrnd(mu1,sigma1,nn);
d2 = mvnrnd(mu2,sigma2,nn);

noise1 = randn(nn,1)*0.1;
noise2 = randn(nn,1)*0.1;

Data=[d1,noise1;d2,noise2];
Labels=[ones(nn,1)*1;ones(nn,1)*2];

figure
gscatter(Data(:,1),Data(:,2),Labels)

figure
gscatter(Data(:,2),Data(:,3),Labels)

figure
gscatter(Data(:,3),Data(:,4),Labels)

figure
gscatter3(Data(:,1),Data(:,3),Data(:,4),Labels)


%%

collabz = {'1';'2';'3';'4'};
testcolumns = [1,2,3,4];
%Train the RBF model
[Fits] = RelativeRBFTrain(Data,Labels,2,'ploton');

%%
testcolumns = [1,2]
collabz = {'1';'2'};
[Variations,BestSep,BestSepClass,Data_All,Diff_All] = RelativeRBFSeperabilityCheck(Data,Labels,testcolumns,collabz,1,'ploton');
BestSep.bestcollabels
BestSep.max
idn = 2;
BestSepClass{idn}.max
BestSepClass{idn}.metric
columnrbf = BestSepClass{idn}.bestcolumns
BestSepClass{idn}.bestcollabels



%% RELIEFF
columnrbf = [1,2];
% 
% figure
% gscatter3(Data(:,columnrbf(1)),Data(:,columnrbf(2)),Data(:,columnrbf(3)),Labels)

DataTest = Data(:,columnrbf);

[RANK,WEIGHT] = relieff(DataTest,Labels,200);

RANK
WEIGHT(RANK)

%% Train RBF discrim
columnrbf = [1,2];
% 
% figure
% gscatter3(Data(:,columnrbf(1)),Data(:,columnrbf(2)),Data(:,columnrbf(3)),Labels)

DataTest = Data(:,columnrbf);


%Train the RBF model
[Fits] = RelativeRBFTrain(DataTest,Labels,3,'ploton');

%try it online
thresh = Fits{end}.ThresholdSep;
[Class,Probability] = RelativeRBFOnline(DataTest,Fits,0.4);

Correct = Class == Labels;
Correct(Class == 0) = [];

accuracy = mean(Correct)