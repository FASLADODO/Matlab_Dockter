function total = signChanges(X)
%signChanges(X); find total number of +/- sign changes in each column of X
%Returns a vector column changes
    X(X<0)=0; X(X>0) = 1;
    temp = diff(X);
    total = sum( abs( temp ));
end