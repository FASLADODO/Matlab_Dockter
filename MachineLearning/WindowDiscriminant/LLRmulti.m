function [LLR,combos] = LLRmulti(Probs)
    [~,SS] = size(Probs);

    combos = nchoosek(1:SS,2);
    [nc,~] = size(combos);
    
    for ii = 1:nc
        LLR(ii,1) = log(sum(Probs(:,combos(ii,1)))/sum(Probs(:,combos(ii,2))));
    end
end