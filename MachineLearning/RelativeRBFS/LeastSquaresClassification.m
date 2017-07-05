function [Class,Accuracy,minError] = LeastSquaresClassification(Data,Labels,LS_Param,modelfunc,outputcolumn)
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]
%modelfunc: gives column order and powers eg:
% modelfunc = @(X) [X(:,1),X(:,2),X(:,1).^2]
% outputcolumn: column for Y = Data(:,outputcolumn)

    [NN,SS] = size(Data);
    cslist = unique(Labels);

    %try classifying with LS
    for cc = 1:length(LS_Param)
        Error(:,cc) = abs(Data(:,outputcolumn) - modelfunc(Data)*LS_Param{cc}.Params);
    end

    %Get class estimate for each data point
    [minError,minidx] = min(Error,[],2);

    %get class for each data point
    Class = cslist(minidx);

    %check accuracy
    Correct = Class == Labels;
    Accuracy = mean(Correct);

end