function [ entropy ] = EntropyCalc( probs )
%For 1D data
    if(sum(probs) ~= 1)
       %disp('Sum of probabilities is not 1'); 
    end
    entropy = 0;
    for kk = 1:length(probs)
        if(probs(kk) ~= 0)
            entropy = entropy - probs(kk)*log2(probs(kk));
        end
    end
end

