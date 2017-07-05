function combos = nchoosekmulti(indices,kmulti)
%does the nchoose k algorithm but with k as an array of values
%k = 1:length(indices) eg
    combos = [];
    for k = kmulti
        combos{k} = nchoosek(indices,k);
    end
end
