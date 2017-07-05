% Searches for argument 'needle' in a 2-D cell array 'haystack'
% returns the column(s) number which contains the needle
% if two return arguments are specified, returns [r c]; the row and column
% number(s) containing the needle
function [varargout] = cellSearch(needle, haystack)

r=[]; c =[];
[m n] = size(haystack);
needle = lower(needle);
for i=1:m
    for j = 1:n
        if ~isempty(strfind(lower(haystack{i, j}), needle))  %ismember(haystack{i, j}, needle) || ...      
            
            r(end+1)=i;
            c(end+1)=j;
        end
    end
end

if (nargout == 2)
    varargout{1}=r;
    varargout{2}=c;
else 
    varargout={c};
end
