function model = RelativeGuassianTrainDiffs(Data,Diff,dataFunc)
%Train parameters and scale to map Data to Difference
%Data: data amtrix with columns as states
%Diff: vector output of ComputeRBFDifference()
%dataFunc: anonymous function to decide on columns
%eg: dataFunc = @(X) [X(:,1), X(:,3)] to use only columns 1 and 3

    %figure out parameters and scale
    dtemp = dataFunc(Data);
    param_diff = pinv(dtemp)*Diff;
    param_diff_scale = max(Diff);
    
    %stash in struct
    model.Params = param_diff;
    model.Scale = param_diff_scale;
    model.Function = dataFunc;
end