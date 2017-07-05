function UnitVector = UnitVector(X)
%computes unit vector for each row of a matrix
%X row wise data matrix, columns are dimensions
[~,SS] = size(X);
N = NormRowWise(X);
UnitVector = X./ repmat(N,1,SS);

end