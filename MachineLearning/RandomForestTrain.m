function [Forest, oobavg] = RandomForestTrain(Data,Labels,cols,ntrees,mfeat)
% random forests implementation
% https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
% Reguires DecisionTreeTrain.m and DecisionTreeOnline.m
% Data: data matrix with rows as samples
% Labels: labels for each sample
% cols: columns labels
% ntrees: number of trees to use
% mfeat: number of random features to use in each tree

trainingratio = 0.63; %two thirds of data

%get size of all data
[NN,MM] = size(Data);
AllIDX = [1:NN]';

%check sizes
if(mfeat >= MM)
   error('mfeat must be lower than the number of features in Data'); 
end

%figure out how many samples to use in each
ndraw = round(trainingratio*NN);
noob = NN - ndraw;

%preallocate
Forest = [];
oobavg = 0;

%all trees
for tt = 1:ntrees
    %figure out which features will get kept
    featcurr = sort(randiunique(MM,mfeat));
    
    %get sample idx
    trainIDX = datasample(AllIDX,ndraw,1,'Replace',true);
    oobIDX = datasample(AllIDX,noob,1,'Replace',true);
    %training data
    trainSample = Data(trainIDX,featcurr);
    trainLabels = Labels(trainIDX,:);
    %out of box data
    oobSample = Data(oobIDX,featcurr);
    oobLabels = Labels(oobIDX,:);
    
    % Make a tree with this subset of the data
    TreeTemp = DecisionTreeTrain(trainSample,trainLabels,cols);
    
    %compute oob error
    LabelEstimate = DecisionTreeOnline(oobSample,TreeTemp);
    wrong = LabelEstimate ~= oobLabels;
    oobError = mean(wrong);
    oobavg = oobavg + oobError;
    
    %Stash it all
    Forest.Tree{tt} = TreeTemp;
    Forest.Tree{tt}.useFeatures = featcurr;
    Forest.Tree{tt}.oobError = oobError;
end
oobavg = oobavg / ntrees; % get average error

end