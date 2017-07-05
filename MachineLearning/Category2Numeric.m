function [class,key] = Category2Numeric(Y)
%convert class labels that are strings to be unqieu integers for math
%requires MapValues() which is a function written by rod

[NN,SS] = size(Y);

if (SS ~=1)
   error('too many columns') 
end

key = unique(Y);
class = MapValues(Y,key,[1:length(key)]);

end