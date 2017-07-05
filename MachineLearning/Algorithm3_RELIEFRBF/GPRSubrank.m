function [Difference,ClassData,DataLabels] = GPRSubrank(Data,Labels)
%just compute difference for each points
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]

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


%loop through all classes
ClassData = [];
Difference = [];
DataLabels = [];

%Find all relative RBFS for class data
for cc = 1:length(cslist)
    %current class and opposite class
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    
    nn = size(DON,1);
    
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    
   
    %compute difference
    Difference = [Difference; ComputeRBFDifference(pwithin,pbetween) ];
    ClassData = [ ClassData; DON];
    DataLabels = [DataLabels; ones(nn,1)*cslist(cc) ];
end

end