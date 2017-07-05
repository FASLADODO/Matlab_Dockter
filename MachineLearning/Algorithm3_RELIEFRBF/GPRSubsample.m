function [SubDifference,SubClassData,SubDataLabels] = GPRSubsample(Data,Labels, ratio)
%sub sample data using RELIEFRBF
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]
%ratio: percentage of points to keep

[NN,SS] = size(Data);

%figure how much to keep
lim = floor(NN*ratio);

%optimal bandwidth (rule of thumb)
sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5));

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

%subset
SubDifference = [];
SubClassData = [];
SubDataLabels = [];

%Find all relative RBFS for class data
for cc = 1:length(cslist)
    %current class and opposite class
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    
    nn = size(DON,1);
    
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    
    lim = floor(nn*ratio);
   
    %compute difference
%     Difference = [Difference; ComputeRBFDifference(pwithin,pbetween) ];
%     ClassData = [ ClassData; DON];
%     DataLabels = [DataLabels; ones(nn,1)*cslist(cc) ];
    Difference = ComputeRBFDifference(pwithin,pbetween) ;
    ClassData = DON;
    DataLabels =  ones(nn,1)*cslist(cc);
    
    % Now sort the differences
    [sortDiff,idd] = sort(Difference,'descend');

    %subsamples
    SubDifference = [SubDifference; sortDiff(1:lim)];
    SubClassData = [SubClassData; ClassData(idd(1:lim),:)];
    SubDataLabels = [SubDataLabels; DataLabels(idd(1:lim),:)];
end

% % Now sort the differences
% [sortDiff,idd] = sort(Difference,'descend');
% 
% %subsamples
% SubDifference = sortDiff(1:lim);
% SubClassData = ClassData(idd(1:lim),:);
% SubDataLabels = DataLabels(idd(1:lim),:);


end