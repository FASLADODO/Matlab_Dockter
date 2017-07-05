function [ClassEstimate,Value] = AdaBoostClassify(Data,Tree)
%Classify using weighted weak learners from adaboost algorithm

    [NN,SS] = size(Data);
    Value = zeros(NN,1);
    %get estimate for each round of boosing, weighted by alpha
    for ii = 1:Tree.Rounds
        Value = Value + AdaBoostDecision(Data,Tree,ii).*Tree.Layer{ii}.Alpha;
    end
    %class esimate is dictated by sign of resultant
    ClassEstimate = sign(Value);
end