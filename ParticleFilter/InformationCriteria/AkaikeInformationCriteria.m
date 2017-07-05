function [AIC,bestcombo] = AkaikeInformationCriteria(Y,X)
%compute akaike information criteria using least squares model
%returns scalar value

[NN,SS] = size(X);

%loop through all possible number of features
for ii = 1:SS

    statepairs = nchoosek(1:SS,ii); %THIS IS THE SHIT!
    [combos, vars] = size(statepairs);
    
    %parameters count
    k = ii + 1;

    for jj = 1:combos
        %get residual sum of squares
        X_sub = X(:,statepairs(jj,:));
        Theta = pinv(X_sub) * Y;
        error = (Y - X_sub*Theta);
        RSS = sum( error.^2 );
        sigma2 = RSS/NN;

        %compute information criteria
        AIC(jj,ii) = 2*k + NN*log(RSS);
        %AIC = -(NN/2)*log(2*pi) - (NN/2)*log(sigma2) - (1/(2*sigma2))*RSS

    end
    
    
    [m, idx] = min(nonzeros(AIC(:,ii)));
    bestcombo{ii} = statepairs(idx,:);
end

    
end