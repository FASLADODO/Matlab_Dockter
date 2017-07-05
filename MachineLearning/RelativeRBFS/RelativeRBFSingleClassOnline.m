function [ClassEstimate,Seperability] = RelativeRBFSingleClassOnline(TrainData,TrainLabels,DataOnline,Direction,Threshold)
%TrainData = class data to train on (rows are samples, columns are dimensions)
%TrainLabels = column of unique class lables with same length as TrainData [1;2...]
%DataOnline = N-D data to try classifying using DataTrain
%Direction = classes to use for probability 
% diff > Threshold = Direction(1)
% diff < Threshold = Direction(1)
% ClassUse = Direction(1)

%class whose forward probability we will use
ClassUse = Direction(1);

%optimal bandwidth (rule of thumb)
sig = norm(std(TrainData));
n = size(TrainData,1);
bw = 1.06*sig*(n^(-1/5));

%get all unique classes
if(isnumeric(TrainLabels))
    cslist = unique(TrainLabels);
    cslookup = strread(num2str(cslist'),'%s');
else
    cslookup = unique(TrainLabels);
    cslist = 1:length(cslookup);
end

%segregate data for each class
%current class and opposite class
DON = TrainData(TrainLabels == ClassUse,:);
DOFF = TrainData(TrainLabels ~= ClassUse,:);

%data sizes
non = size(DON,1);
noff = size(DOFF,1);

%get within and between class
within_us = rbfpdist2(DataOnline,DON,bw);
between_us = rbfpdist2(DataOnline,DOFF,bw);

%scale it
pwithin = double(within_us/non);
pbetween = double(between_us/noff);
    
%compute difference
Seperability = ComputeRBFDifference(pwithin,pbetween);

%Estimate class
ClassEstimate = Direction( (Seperability < Threshold) + 1);

end