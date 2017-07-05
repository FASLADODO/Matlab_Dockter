function [Difference] = RelativeRBFSingleClass(Data,Labels,ClassUse)
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]
%ClassUse = class to use for probability

[NN,SS] = size(Data);

%optimal bandwidth (rule of thumb)
sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))

%get all unique classes
if(isnumeric(Labels))
    cslist = unique(Labels);
    cslookup = strread(num2str(cslist'),'%s');
else
    cslookup = unique(Labels);
    cslist = 1:length(cslookup);
end

%segregate data for each class
%current class and opposite class
DON = Data(Labels == ClassUse,:);
DOFF = Data(Labels ~= ClassUse,:);

%data sizes
non = size(DON,1);
noff = size(DOFF,1);

%get within and between class
within_us = rbfpdist2(Data,DON,bw);
between_us = rbfpdist2(Data,DOFF,bw);

%scale it
pwithin = double(within_us/non);
pbetween = double(between_us/noff);
    
%compute difference
Difference = ComputeRBFDifference(pwithin,pbetween);

end