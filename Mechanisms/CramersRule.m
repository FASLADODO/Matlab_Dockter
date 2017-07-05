function [vals] = CramersRule(A,B)
    vals = [];
    [rows, cols] = size(A);
    
    for i = 1:rows
        tempA = A;
        tempA(:,i) = B;
        tempval = det(tempA)/det(A);
        vals = [vals, tempval];
    end


end