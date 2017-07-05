function N = NormRowWise(X)
%Compute 2-norm for each row of the matrix X

N = sqrt(sum(abs(X).^2,2));

end