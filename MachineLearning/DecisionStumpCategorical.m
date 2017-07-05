function [Thresh,Direction,IG] = DecisionStumpCategorical(X,Labels)
%Computes optimal category for X numerical data to classify with
%X is a single column of data categorical
%W are the current weights for each samples
%Labels are class labels for each sample
%Returns:
%Thresh: optimal class
%best score is the sum of weights of incorrectly classified data
%IG: information gain from classify with this best threshold
%http://www.autonlab.org/_media/tutorials/dtree18.pdf

[NN,SS] = size(X);
if(SS > 1)
   warning('X has too many columns, using only the first')
   X = X(:,1);
end

%groups
groupings = unique(X);
cslist = unique(Labels);

%initialize it
Thresh = 0;
IG = 0;
Direction = cslist';

%loop thorugh all possible thresholds
for ii = 1:length(groupings)
    tempdir = cslist;
    
    %try a grouping
    isabove = X == groupings(ii);
    
    %see which class that is
    ClassEst = tempdir(isabove+1);
    
    %see whats incorrect above and below threshold
    wrong = ClassEst ~= Labels;
    
    %see if we have to flip ie: the hypothesis is backwards
    if(mean(wrong) > 0.5)
        tempdir = flipud(tempdir);
        ClassEst = tempdir(isabove+1);
    end
    
    %compute information gain
    IG_Temp = InformationGain(Labels,ClassEst);
    
    %should we keep it
    if(IG_Temp > IG)
       Thresh = groupings(ii);
       IG = IG_Temp;
       Direction = tempdir';
    end
end
