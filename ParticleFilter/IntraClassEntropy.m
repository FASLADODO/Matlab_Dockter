function [ entropy ] = IntraClassEntropy( class1, class2 )
    combined = (class1 + class2)/2;
    entropy = 0;
    for kk = 1:length(combined)
        if(combined(kk) ~= 0)
            entropy = entropy - combined(kk)*log2(combined(kk));
        end
    end

end

