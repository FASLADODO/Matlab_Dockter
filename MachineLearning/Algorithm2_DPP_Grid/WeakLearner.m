function [Thresh,Acc,Dir] = WeakLearner(X,Labels,cslist)
%Computes optimal threshold for X numerical data to classify with
%X is a single column of data
%Labels are class labels for each sample
%Returns:
%Thresh: optimal threshold
%Acc: is the sum of weights of incorrectly classified data
%IG: information gain from classify with this best threshold
%http://www.autonlab.org/_media/tutorials/dtree18.pdf

[NN,SS] = size(X);
if(SS > 1)
   warning('X has too many columns, using only the first')
   X = X(:,1);
end

%Fill in  sorted data
stepratio = 0.5;

%sort em
XSort = sort(X);

%steps
DX = diff(XSort);

%initialize it
Thresh = min(X);
acc_Temp = 0;
Acc = 0;
Dir = cslist';

%loop thorugh all possible thresholds
for ii = 1:NN-1
    %try a threshold
    threshtemp = XSort(ii,:) + (DX(ii,:)*stepratio);
    isbelow = X < threshtemp;

    %classify based on this threshold
    ClassEst1 = cslist(isbelow + 1);
    ClassEst2 = cslist(~isbelow + 1);
    
    %see how many are correct
    Corr1 = ClassEst1 == Labels;
    Corr2 = ClassEst2 == Labels;
    
    %accuracies
    Acc1 = mean(Corr1);
    Acc2 = mean(Corr2);
    
    %figure out which direction gave which class
    [acc_Temp,id] = max([Acc1,Acc2]);
    
    %should we keep it
    if(acc_Temp > Acc)
       Thresh = threshtemp;
       Acc = acc_Temp;
       if(id == 1)
           Dir = cslist';
       else
           Dir = fliplr(cslist');
       end
       
    end
end


end