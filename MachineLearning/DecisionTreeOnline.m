function [ClassEstimate] = DecisionTreeOnline(X,Tree)
%X is data matrix with rows as a sample
%Tree is a decision tree trained by DecisionTreeTrain

[NN,SS] = size(X);

for ii = 1:NN
    nind = 1; %always start from the top of the tree
    prevind = 0;
    X_Current = X(ii,:);
    while(1)
        if(~isempty(Tree.Leaves{nind}.Resultant) )
            ClassEstimate(ii,:) = Tree.Leaves{nind}.Resultant;
            break;
        end
        %make sure we aren't stuck in an infinite loop
        if(nind == prevind)
            ClassEstimate(ii,:) = [];
            break;
        end
        prevind = nind;
        
        %get tree parameters at this node
        feat = Tree.Leaves{nind}.Feature;
        thresh = Tree.Leaves{nind}.Threshold;
        
        %get the current feature
        X_Feature = X_Current(feat);
        %see where we go next
        if(X_Feature < thresh)
            nind = Tree.Leaves{nind}.ltindex;
        elseif(X_Feature >= thresh)
            nind = Tree.Leaves{nind}.gtindex;
        end
    end
end


end