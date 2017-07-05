function mask = inRange(Data,minD,maxD)
%check which rows are in range given limits
%Data: rows are samples,column are dimensions
%min = min(Data,[],1)
%max = max(Data,[],1)

%EG:
% Data = rand(10,2)
% minD = min(Data,[],1)
% maxD = max(Data,[],1)
% mask =  inRange(Data,minD,maxD)

[NN,~] = size(Data);

%Check if each row is between min and max in all dimensions
AllColumns = [(Data >= repmat(minD,NN,1)), (Data <= repmat(maxD,NN,1))];

%check for outside in any
mask = logical(min(AllColumns,[],2));

end