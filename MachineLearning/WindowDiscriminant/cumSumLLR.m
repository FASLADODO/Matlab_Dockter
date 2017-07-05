function [LLR,combos] = cumSumLLR(Probs)

    [~,SS] = size(Probs);

    combos = nchoosek(1:SS,2);
    [nc,~] = size(combos);
    
    for ii = 1:nc
        LLR(:,ii) = log(cumsum(Probs(:,combos(ii,1))) ./ cumsum(Probs(:,combos(ii,2))));
    end
end