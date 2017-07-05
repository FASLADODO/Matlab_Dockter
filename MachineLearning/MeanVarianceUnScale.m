function [Xhat] = MeanVarianceUnScale(X,shift,scale)
%scale each dimension of X to be zero mean and var = 1
%X = rows as samples, columns as dimensions

%get size
[NN,~] = size(X);

%scale and shift
Xhat = X + repmat(shift,NN,1);
Xhat = Xhat .* repmat(scale,NN,1);

end