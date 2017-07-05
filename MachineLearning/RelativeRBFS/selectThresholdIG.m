function [Thresh,IG,AllThresh] = selectThresholdIG(Difference,Labels,Direction)
%Difference is a column vector of seperabilities or probabilities for each
%This comes from RelativeRBFSingleClass()
%Labels are class membership
%Direction is a nx1 vector is [1,2] where data > thresh is class 1,
%data < thresh is class 2.
%Direction(1) = class whose probabilities are used ie:
%classon = Difference(Labels == Direction(1),:) (should be central data)

classon = Direction(1);

%get class specific seperabilities and sort them
DiffClass = sort(Difference(Labels == classon,:));
DX = diff(DiffClass);


IG = 0;
Thresh = 0;
AllThresh = [];

%loop through all difference data for class on (test thresholds)
for ii = 1:length(DiffClass)-1 
    tempthresh = DiffClass(ii) + (DX(ii) / 2);
    
    %see if we're above that threshold
    isabove = Difference < tempthresh;
    
    %estimate class basedo n is above (1 is above, 2 is below)
    classest = Direction(isabove+1);
    
    %compute information gain given classification using current threshold
    IG_Temp = InformationGain(Labels,classest);
    
    %store IG and thresh
    AllThresh = [AllThresh; tempthresh, IG_Temp];
    
    if(IG_Temp > IG)
       IG = IG_Temp;
       Thresh = tempthresh;
    end
end


end