function [Model] = RelativeRBFWeightedLS(Data,Labels,modelfunc,outputcolumn)
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


%Find all relative RBFS for class data
for cc = 1:length(cslist)
    %current class and opposite class
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    
    nn = size(DON,1);
    
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    
    %get weights from difference
    Difference = ComputeRBFDifference(pwithin,pbetween);
    
    W = diag(Difference);
    
    %get data columns
    DTrain = modelfunc(DON);
    %get output column
    YTrain = DON(:,outputcolumn);
    %weighted least squares
    Params = inv(DTrain'*W*DTrain)*DTrain'*W*YTrain;
    
    %
    Model{cc}.Weights = Difference;
    Model{cc}.Params = Params;
    Model{cc}.Class = cslist(cc);
end