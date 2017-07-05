function XPRIME = NormalizeFeatures(X)
%normalizes each feature using mean and variance
%results should scale to -1-1 for all features
%X is a data matrix, rows are samples, columns are features
%Returns XPRIME which is the same size as X

%get size
[NN,SS] = size(X);
XPRIME = zeros(NN,SS);

%loop through all features
for ii = 1:SS
   dtemp = X(:,ii); %current dimensions data
   mu = (1/NN)*sum(dtemp); %mean
   sigma = (1/NN)*sum((dtemp - mu).^2); %variance
   XPRIME(:,ii) = (1/sigma).*(dtemp - mu); % scale it
end

end
