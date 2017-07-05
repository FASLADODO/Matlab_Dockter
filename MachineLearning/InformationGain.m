function IG = InformationGain(Y,X,threshold)
%compute the estimated information gain Y given X
%"how much entropy is saved given knowledge of X"
%X can be numerical or categorical 
%if X is numerical supply threshold to get 0/1 true false
%if X is categorical don't supply threshold
%threshold optional scalar
%slide 11
%http://www.cs.cmu.edu/~guestrin/Class/10701-S06/Handouts/recitations/recitation-decision_trees-adaboost-02-09-2006.ppt
    
%First check if Y is string
if(iscell(Y))
    [class,~] = Category2Numeric(Y);
    Y = class;
end
cslist = unique(Y);

%then get entropy of Y
H_Y = EntropySample(Y);

%Next check if X needs a threshold
if(nargin == 3)
    %get conditional entropies using threshold to group X
    H_Y_X = ConditionalEntropy(Y,X,threshold);
else
    %get conditional entropies for categorical X
    H_Y_X = ConditionalEntropy(Y,X);
end

%information gain
IG = H_Y - H_Y_X;


end






