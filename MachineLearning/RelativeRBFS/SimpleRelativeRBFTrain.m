function [Difference,ClassData,Model] = SimpleRelativeRBFTrain(Data,Labels)
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]
%kmeanlimit = estimate of max number of clusters for each class
%ploton = 'ploton' or 'plotoff'

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

%Find all relative RBFS for class data
for cc = 1:length(cslist)
    %current class and opposite class
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    
    nn = size(DON,1);
    
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    
   
    %compute difference
    Difference{cc} = ComputeRBFDifference(pwithin,pbetween);
    ClassData{cc} = DON;
end

%loop through all classes
Model = [];
Model{1}.cslist = cslist;

%Find the good data
for cc = 1:length(Difference)
    %current class
    DON = ClassData{cc};
    nn = size(DON,1);
    if(nn > 0)
        Diff = Difference{cc};
        %thresh = 0.1*max(Diff);
        %thresh = 0.2*prctile(Diff,95);
        thresh = 0.1*prctile(Diff,95);

        %find all data above a threshold
        idx = find(Diff > thresh);

        %stash it
        GoodData = DON(idx,:);
        
        %retry difference with good data
        pwithin = rbfpdist2(GoodData,GoodData,bw) / nn;
        threshgood = prctile(pwithin,5);
        
        %Make a model
        Model{cc}.gooddata = GoodData;
        Model{cc}.prob = pwithin;
        Model{cc}.thresh = threshgood;
    end
end

end