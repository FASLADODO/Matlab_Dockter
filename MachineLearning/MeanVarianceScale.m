function [Xhat,shift,scale] = MeanVarianceScale(X,shift,scale)
%scale each dimension of X to be zero mean and var = 1
%X = rows as samples, columns as dimensions

%get size
[NN,~] = size(X);

if(nargin == 1)
    %get mean and variances
    shift = mean(X);
    scale = std(X);
end

%scale and shift
Xhat = X - repmat(shift,NN,1);
Xhat = Xhat ./ repmat(scale,NN,1);

end