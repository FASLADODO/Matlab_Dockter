function Model = LLSRtrain(X,Y,Labels)
%Train a least squares model with Logistic regression for confidence
%measure
%X: data matrix with columns as dimensions
%Y: Output value of X when mapped through model
%Labels: unique label for each class in training

%Returns Model which holds least squares params and logistic regression
%parameters
%Model.LSParams: cell array with least squares parameters
%Model.LogReg: parameters for logistic regression

    classes = unique(Labels);

    [NN,SS] = size(X);
    [NY,SY] = size(Y);
    
    if(NN ~= NY)
       error('X and Y must be same length')
    end

    for cc = 1:length(classes)
        dtemp = X(Labels == classes(cc) , :);
        ytemp = Y(Labels == classes(cc) , :);

        phitemp = pinv(dtemp)*ytemp;

        Model.LSParams{cc} = phitemp;

    end
    
    %concanenate data
    LogData = [X,Y];
    
    %Now get logistic parameters
    Model.LogReg =  LogRegTrain(LogData,Labels);
    
    %run logistic regression through
    [PALL,~] = LogRegOnline(LogData,Model.LogReg );
    [LLR,~] = cumSumLLR(PALL);
    
    Model.LogScale = max(abs(LLR));
end