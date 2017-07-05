function [BIC_Results] = BICRMSE(Data,X_World,Y_World)
%checks all combinations of all variables in Data and sees whcih
%combinations give the lowest error

    [NN,SS] = size(Data);

    BIC_Results = [];

    for dd = 1:SS
        %all possible combinations with this number of dimensions
        combinations = nchoosek([1:SS],dd);

        %best rmse
        best_rmse = inf;
        best_columns = [];

        %loop through all n choose k
        for ii = 1:size(combinations,1)
            %columns to test
            columns_use = combinations(ii,:);

            %get this data
            DataTest = Data(:,columns_use);

            %least squares params
            paramsX = pinv(DataTest)*X_World;
            paramsY = pinv(DataTest)*Y_World;

            %test it
            X_Estimate = DataTest*paramsX;
            Y_Estimate = DataTest*paramsY;

            %compute error and rmse
            Error = norm([X_World-X_Estimate,Y_World-Y_Estimate]);
            rmse = sqrt(mean(Error));

            %stash it
            BIC_Results.Dimensions{dd}.AllVariationsColumns(ii,:) = columns_use;
            BIC_Results.Dimensions{dd}.AllVariationsRMSE(ii,:)= rmse;

            %find best 
            if(rmse < best_rmse)
                best_rmse = rmse;
                best_columns = columns_use;
            end
        end
        %stash it
        BIC_Results.Dimensions{dd}.BestVariation.columns = best_columns;
        BIC_Results.Dimensions{dd}.BestVariation.rmse = best_rmse;
        
        %for easy plotting
        BIC_Results.BestRMSE(dd,:) = best_rmse;
    end

end
