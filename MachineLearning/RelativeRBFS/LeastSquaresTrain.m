function [LS_Param] = LeastSquaresTrain(Data,Labels,modelfunc,outputcolumn)
%Data = class data to train on (rows are samples, columns are dimensions)
%Labels = column of unique class lables with same length as Data [1;2...]
%modelfunc: gives column order and powers eg:
% modelfunc = @(X) [X(:,1),X(:,2),X(:,1).^2]
% outputcolumn: column for Y = Data(:,outputcolumn)

    [NN,SS] = size(Data);
    cslist = unique(Labels);

    %find least squares params
    for cc = 1:length(cslist)
        Dtemp = Data(Labels == cslist(cc),:);
        %Get training info
        DTrain = modelfunc(Dtemp);
        YTrain = Dtemp(:,outputcolumn);

        %calulcate least squares params
        LS_Param{cc}.Params = pinv(DTrain)*YTrain;
        
        LS_Param{cc}.Class = cslist(cc);
    end
    
end
