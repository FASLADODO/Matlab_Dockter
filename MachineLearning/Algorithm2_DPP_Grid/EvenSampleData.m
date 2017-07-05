function [DataEven,LabelsEven] = EvenSampleData(Data,Labels)
%give equal samples to all classes

% number of classes and data sizes
cslist = unique(Labels);
[NN,SS] = size(Data);

%figure out samples from each class
for cc = 1:length(cslist)
    ltemp(cc) = sum(Labels == cslist(cc));
end

%which to subsample
subsize = min(ltemp);

DataEven = [];
LabelsEven = [];
%Subsample it
for cc = 1:length(cslist)
    dtemp = Data(Labels == cslist(cc),:);
    dnew = datasample(dtemp,subsize);
    DataEven = [DataEven; dnew];
    LabelsEven = [LabelsEven; ones(subsize,1)*cslist(cc) ];
end


end