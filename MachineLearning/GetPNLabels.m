function [X,Labels] = GetPNLabels(X,Labels,keep)
%if you want to make your labels [-1,1]
%and you don't care if you lose classes
%keep = [1,2] will keep unique label 1 and 2

[Labels,key] = Category2Numeric(Labels)

%get [-1,1] vals
for ii = 1:length(key)
    idx = find(Labels == ii);
    if(ii == keep(1))
        Labels(idx) =  -1;
    elseif(ii == keep(2))
       Labels(idx) =  1; 
    else
        Labels(idx) =  []; 
        X(idx,:) =  []; 
    end
end