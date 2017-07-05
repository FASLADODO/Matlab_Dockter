function [eigenvals] = CCAbasic(X,Labels)
%Canonical Correlation coefficient
%https://en.wikipedia.org/wiki/Canonical_correlation

%TEST CCA

%X = Data matrix with columns as states, rows as observations
%Labels = class labels
%eigenvals = eigenvalues


[NN,SS] = size(X);
if(SS > NN)
    X = X'; 
    [NN,SS] = size(X);
end


cslist = unique(Labels);

for ii = 1:length(cslist)
    idx = Labels == cslist(ii);
    X_Class{ii} = X(idx,:);
end

% make it zero mean data
X_Bar = X - repmat(mean(X),NN,1); 

%covariance of mean adjusted data
covX = cov(X_Bar);

%get eigenvec and eigenvals
[eigV,eigD] = eig(covX);

%Sort the eigenvalues
[eigD_sort,idxV] = sort(diag(eigD)','descend');

%sort eigenvectors according to eigenvalues
eigV_sort = eigV(:,idxV);
eigV_sort

%Get feature vector for the number of components we wish to keep
FeatureVector = eigV_sort(:,[1:NC])';

%now transform to get Data in PCA coordinate system
X_PCA = FeatureVector*X_Bar';
X_PCA = X_PCA'; %columns are easier


end