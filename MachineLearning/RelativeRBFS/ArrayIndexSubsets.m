function arrayIndexSets = ArrayIndexSubsets(totalIdx,limit)
%will return a matrix of indices that loops the total number of indices
%into subsets 1:limit ie:
%totalIdx = 8;
%limit = 3;
%arrayIndexSets = [1,2,3;
%                 4,5,6;
%                 6,7,8];


%get divisor
divisor = totalIdx / limit;

%get rounding
rounded = floor(divisor);

%remainder
remain = divisor - rounded;

%initialize
arrayIndexSets = []; 
prevstart = 0;

%loop through the number of rounded ones
for ii = 1:rounded
    temp = 1:limit; %indices
    newidx = temp + prevstart; %increment
    arrayIndexSets(ii,:) = newidx; %store in this row
    prevstart = max(newidx);
end

%store the last row
if(remain > 0)
    temp = (totalIdx - limit + 1):totalIdx;
    arrayIndexSets(end+1,:) = temp;
end

end

