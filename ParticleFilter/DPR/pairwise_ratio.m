function class = pairwise_ratio(probs,thresh)
    for ii = 1:length(probs)
       sumden = 0;
       %sum all other vals
       for jj = 1:length(probs)
           if(ii ~= jj)
              sumden = sumden + probs(jj); 
           end
       end
       %compute ratios
       if(sumden == 0)
        ratio(ii) = probs(ii);
       else
        ratio(ii) = probs(ii) / sumden;
       end
    end
    %get max ratio
    [val,idx] = max(ratio);
    if(idx > thresh)
        class = idx;
    else
        class = 0;
    end
end