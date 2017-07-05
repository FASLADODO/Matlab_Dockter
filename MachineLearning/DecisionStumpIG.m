function [Thresh,IG] = DecisionStumpIG(X,Labels)
%Computes optimal threshold for X data using information gain
%X is a single column of data
%Labels are class labels for each sample
%Returns:
%Thresh: optimal threshold
%IG: information gain from classify with the best threshold
%http://www.autonlab.org/_media/tutorials/dtree18.pdf

[NN,SS] = size(X);
if(SS > 1)
   warning('X has too many columns, using only the first')
   X = X(:,1);
end

%Fill in missing info if minimal function is called
cslist = unique(Labels);

%figure out which class has lower mean
for ii = 1:length(cslist)
    dtemp = X(Labels == cslist(ii),:);
    mu(ii) = mean(dtemp);
end
%figure out which class has the lowest mean
[~,muid] = sort(mu);
Direction = cslist(muid);

%sort em
XSort = sort(X);

%steps
DX = diff(XSort);

%step size default is 2 ie X + DX/2
stepratio = 2;

%initialize it
Thresh = 0;
IG = 0;

%loop thorugh all possible thresholds
for ii = 1:NN-1
    %try a threshold
    threshtemp = XSort(ii,:) + (DX(ii,:)/stepratio);
    isabove = X >= threshtemp;
    
    %see which class that is
    ClassEst = Direction(isabove+1);
    
    %compute information gain
    IG_Temp = InformationGain(Labels,ClassEst);
    
    %should we keep it
    if(IG_Temp > IG)
       Thresh = threshtemp;
       IG = IG_Temp;
    end
end


end