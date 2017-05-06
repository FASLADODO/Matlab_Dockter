function [X_Sort,IDX_Sort] = InsertionSort(X)

    %initialize
    X_Sort = [];
    IDX_Sort = [];
    
    [NN,SS] = size(X); %vector size
    if(SS > 1)
       error('input must be single column vector') 
    end
    
    %loop through original vecotr
    for ii = 1:NN

        ins = 0;
        %loop through new sorted vector
        for jj = length(X_Sort):-1:1
            if(X(ii) >= X_Sort(jj))
                %if it fits somewhere, put in the array
                X_Sort = [X_Sort(1:jj),X(ii),X_Sort(jj+1:end)];
                IDX_Sort = [IDX_Sort(1:jj),ii,IDX_Sort(jj+1:end)];
                ins = 1;
                break;
            end
        end
        
        if(ins == 0)
           %just tack it at beggining
           X_Sort = [X(ii),X_Sort];
           IDX_Sort = [ii,IDX_Sort];
        end
    end
    
    X_Sort = X_Sort';
    IDX_Sort = IDX_Sort';
end