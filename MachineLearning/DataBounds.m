function Bounds = DataBounds(X)
%returns 2 rows and N columns with the upper and lower bounds for each
%dimension of matrix X
    Bounds = [min(X); max(X)];

end