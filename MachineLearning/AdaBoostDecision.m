function ClassEstimate = AdaBoostDecision(Data,Tree,idt)
%given a tree and some data compute the class estimate
%only for single tree classification
%for use with adaboost

    X = Data(:,Tree.Layer{idt}.Dimension);
    Thresh = Tree.Layer{idt}.Threshold;
    Direction = Tree.Layer{idt}.Direction;
    %using decision stumps
    ClassEstimate = DecisionStumpOnline(X,Thresh,Direction);
end