function [ClassEst, ClassEstCum, Confidence,ScoreAll] = DLS_Online(X, Y, parameters)
% Estimate the true class based on DLS parameters

[NN,SS] = size(X);
[~,classes] = size(parameters);

%Loop through all classes, get errors
for cc = 1:classes
    %compute error for each parameter set
    errtemp = abs( Y - (X * parameters(:,cc) ) );
    delta(:,cc) = errtemp;
    sumcum(:,cc) = cumsum(errtemp,1);
end

ScoreAll = [];
Confidence = [];
%Get classes from errors
for cc = 1:classes
    num_ratio = sumcum(:,cc);
    den_ratio = max(sumcum,[],2);

    %loop through all pairwise things
    for jj = 1:classes
        if cc ~= jj
            num_ratio = num_ratio - sumcum(:,jj);
            %den_ratio = den_ratio + sumcum(:,jj);
        end
    end
    num_ratio = abs(num_ratio);

    %compute ratios
    if(den_ratio ~= 0 )
        w_ratio = num_ratio./den_ratio;
    else
        den_ratio = den_ratio + eps;
        w_ratio = num_ratio./den_ratio;
    end

    Confidence(:,cc) = w_ratio;
    
end

[bestper,~] = max(delta,[],2);
[worstper,~] = min(delta,[],2);
ScoreAll = abs(bestper-worstper) ./ bestper;

%get overall class
[~,ClassEstCum] = min(sumcum,[],2);
[~,ClassEst] = min(sum(delta));


end