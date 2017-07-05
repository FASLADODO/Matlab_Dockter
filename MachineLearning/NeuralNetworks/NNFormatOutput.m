function [YBinary,cslist,mapping] = NNFormatOutput(Y)
%Format class labels to be 0 or 1
%creates extra columns for more than two classes

    %sort unique values into columns of 0s and 1s
    cslist = unique(Y);
    if(length(cslist) == 2)
        %just two classes
        mapping = [0;1];
        YBinary = zeros(length(Y),1);
        for ii = 1:length(cslist)
            YBinary(Y == cslist(ii),:) = mapping(ii);
        end
    else
        %more than two classes
        mapping = dec2bin(1:length(cslist)) - '0';
        clss = size(mapping,2);
        YBinary = zeros(length(Y),clss);
        for ii = 1:length(cslist)
            if(iscell(Y)) %if classes are strings

                %YBinary(strcmp([Y], cslist(ii)),:) = mapping(ii,:);
            else %if classes are just numbers
                idx = find(Y == cslist(ii));
                YBinary(idx,:) = repmat(mapping(ii,:),length(idx),1);
            end
        end
    end

end