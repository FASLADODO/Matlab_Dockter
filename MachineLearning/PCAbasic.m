function [coeff,score,latent,X_PCA] = PCAbasic(X,W)
%http://www.cs.otago.ac.nz/cosc453/student_tutorials/principal_components.pdf

%X = Data matrix with columns as states, rows as observations
%W = weights of each dimension ie w = 1./var(X);
%coeff = returns the eigenvector weight matrix to get data new PCA coordinates
% score is just the data scaled by the weights
% latent is the eigen values
% X_PCA is the mapped data

%NOw start actual PCA
[NN,SS] = size(X);
if(SS > NN)
    X = X'; 
    [NN,SS] = size(X);
end

% make it zero mean data
mu = mean(X);
X = bsxfun(@minus,X,mu);

%see if we need to scale
if(nargin == 2)
    %weight each feature
    X = bsxfun(@times,X,sqrt(W));
end

%covariance of mean adjusted data
covX = (X'*X)/(NN-1);

%get eigenvec and eigenvals
[coeff,eigvals] = eig(covX);

%Sort the eigenvalues
[eigValues,idx] = sort(diag(eigvals)','descend');

%sort eigenvectors according to eigenvalues
coeff = coeff(:,idx);

%Get feature vector scaled by weights
if(nargin == 2)
    coeff = bsxfun(@times, coeff, 1./sqrt(W)');
end

%now get scores
score = X/coeff';
latent = eigValues;

%transform
X_PCA = X' * coeff;


end