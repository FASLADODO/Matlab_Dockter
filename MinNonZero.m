function [M,IDX] = MinNonZero(X,dim)
    [M,IDX] = min(X(X > 0),[],dim);
end