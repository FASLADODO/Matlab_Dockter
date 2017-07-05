function [rank,weights] = relieffBasic(X,Labels,K)
%implement simple version of relieff algorithm
%https://en.wikipedia.org/wiki/Relief_(feature_selection)
%DOES NOT WORK AS GOOD AS BUILTIN
    
    %have to scale each state evenly
    % Find max and min for every state
    Xmax = max(X);
    Xmin = min(X);
    Xdiff = abs(Xmax-Xmin);
    
    % Exclude single-valued attributes
    isOneValue = Xdiff < eps(Xmax);
    %X(:,isOneValue) = [];
    %Xdiff(isOneValue) = [];
    
    %get data size
    [NN,SS ] = size(X);
    cslist = unique(Labels);
    classmax = sum(cslist);

    % Scale and center the states
    X = bsxfun(@rdivide,bsxfun(@minus,X,mean(X)),Xdiff);
    
    classPrb = [];
    %get data for each class
    for ii = 1:length(cslist)
        idx = Labels == cslist(ii);
        Datat = X(idx,:);
        [Nt{ii},~] = size(Datat);
        if(K > Nt{ii})
            warning('K is too large')
            K = Nt{ii};
        end
        classPrb(ii) = Nt{ii}/NN;
    end
    
    
    weights = zeros(1,SS);
    %sort through all classes
    for ii = 1:length(cslist)
        numupdates = round(Nt{ii}/2);
        rndIdx = randsample(Nt{ii},numupdates);
        %sort through all features
        for kk = 1:SS
           if(~isOneValue(kk))
               idxon = Labels == cslist(ii);
               idxoff = Labels ~= cslist(ii);
               DataOn = X(idxon,:);
               DataOff = X(idxoff,:);
               %get distances to all hits and all misses
               temphit = pdist2(DataOn(rndIdx,kk),DataOn(rndIdx,kk) );
               tempmiss = pdist2(DataOn(rndIdx,kk),DataOff(:,kk) );
               %find nearest neighbors in actual class and opposite class
               sorthit = sort(temphit,2);
               sortmiss = sort(tempmiss,2);
               %average the distances to nearest neighbors
               nearhit =  mean( mean(sorthit(:,2:1+K),2 ) ) ; % *classPrb(ii);
               nearmiss = mean( mean(sortmiss(:,1:K),2  ) ) / (1-classPrb(ii));
               weights(kk) = weights(kk)-nearhit + nearmiss;
           end
        end
    end

    [ws,rank] = sort(weights,'descend');
end