function [ClassEstimate] = RandomForestOnline(X,Forest)
% random forests classification
% https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
%X is data matrix with rows as a sample
%Forest is a decision tree ensemble trained by RandomForestTrain

[NN,~] = size(X);
ClassEstimate = zeros(NN,1);

for ii = 1:NN
    %grab current sample
    X_Current = X(ii,:);
    classEstTrees = [];
    for tt = 1:length(Forest.Tree)
        %get the current tree
        Tree = Forest.Tree{tt};
        %reset indices for tree
        nind = 1; %always start from the top of the tree
        prevind = 0;
        while(1)
            if(~isempty(Tree.Leaves{nind}.Resultant) )
                classEstTrees(tt,:) = Tree.Leaves{nind}.Resultant;
                break;
            end
            %make sure we aren't stuck in an infinite loop
            if(nind == prevind)
                classEstTrees(tt,:) = [];
                break;
            end
            prevind = nind;

            %get tree parameters at this node
            feat = Tree.useFeatures(Tree.Leaves{nind}.Feature);
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
    %now try classifying using resultant class
    ClassEstimate(ii,:) = mode(classEstTrees);
end


end