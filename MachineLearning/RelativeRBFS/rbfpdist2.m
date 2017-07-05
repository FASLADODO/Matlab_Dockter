function p = rbfpdist2(Data1,Data2,bw)
    %get rbf at each point in data1 given data2 using a for loop to avoid
    %huge square matrices
    
    usefor = true;
    
    if(usefor)
        %get probability at each point in grid
        %same as doing each row of square matrix seperately
        parfor ii = 1:size(Data1,1)
            dist = pdist2(Data1(ii,:),Data2);
            r = exp( -(dist.^2)./(2*bw^2) );
            p(ii,:) = sum(r,2);
        end
    else
        %get probability in big old square distance matrix
        %(fails for big data sets)
        dist = pdist2(Data1,Data2);
        r = exp( -(dist.^2)./(2*bw^2) );
        p = sum(r,2);
    end
    
    %deal with zeros so we don't divide by them
    p(p < 2*eps) = 2*eps;
end