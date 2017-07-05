function P = correlationPearson(X,Y)
%compute correlation using pearsons for variables (column of X and Y)
%https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient

[NN,SS] = size(X);
[NY,SY] = size(Y);

if(SY ~= SS)
   error('X and Y must have the same amount of columns') 
end

if(NY ~= NN)
   error('X and Y must have the same amount of rows') 
end

if (NN < 10)
   error('need moar data') 
end

% create square correlation matrix
for ii = 1:SS
    for jj = 1:SS
        %get state data
        xt = X(:,ii); 
        yt = Y(:,jj); 
        
        %zero mean shift
        Ex = xt - mean(xt);
        Ey = yt - mean(yt);

        %cov(x,y) / var(x)*var(y)
        P(ii,jj) = (sum(Ex.*Ey)) / (sqrt(sum(Ex.^2))*sqrt(sum(Ey.^2)));
    end
end



end
