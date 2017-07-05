function [Y] = NNUnformatOutput(YBinary,cslist,mapping)
%Format class labels back to original
%creates extra columns for more than two classes

    [NNo,~] = size(YBinary);
    [NNm,~] = size(YBinary);
    Y = zeros(NNo,1);
    for tt = 1:NNo
       temp =  YBinary(tt,:);
       for mm = 1:NNm
           if(temp == mapping(mm,:))
               Y(tt,:) = cslist(mm);
               break;
           end
       end
    end
end