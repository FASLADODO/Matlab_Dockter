function [bestModel, minRMS, bestParams, bestCoords, AllErr] = ModelSelection(X,Y)
    %find the best Least Squares model for data X given output Y
    %X should ideally have every columns that were may not be relevent to Y

    %seperate data
    [idxtrain, idxval, idxtest]  = dividerand(size(X,1),0.7,0.15,0.15);
    
    %return sub sampled vector
    trainX = X(idxtrain,:);
    valX= X(idxval,:);
    testX = X(idxtest,:);
    trainY = Y(idxtrain,:);
    valY= Y(idxval,:);
    testY = Y(idxtest,:);
    


    %Get sizes of data sets
    [NNm,SSm] = size(trainX);
    [NNt,SSt] = size(testX);
    
    %get indices of models
    Models = 1:SSm;

    %for storage
    minRMS = inf;
    bestModel = [];
    bestCoords = [];
    bestParams = [];

    %Loop through each possible number of paramaeters and each model
    AllErr = [];
    for mm = 1:length(Models)
        %get all possible variations for models with that number of terms
        vars = nchoosek(Models,Models(mm));

        [rv,cv] =size(vars);
        for vv = 1:rv
            combv = vars(vv,:);

            %train parameters
            Xon = trainX(:,[combv] );
            params = pinv(Xon)*trainY;

            %get test data set
            Xtest = testX(:,[combv] );
            YBar = Xtest*params;

            %compute rms error for that model
            Errtemp = (1/NNt) * sum( abs(YBar - testY) .^2 );

            %store errors and models
            AllErr(mm,vv) = Errtemp;
            Modelstore{mm,vv} = combv;

            %check if we ofund the best yet
            if(Errtemp < minRMS)
                minRMS = Errtemp;
                bestModel = combv;
                bestCoords = [mm,vv];
                bestParams = params;
            end
        end

    end

end