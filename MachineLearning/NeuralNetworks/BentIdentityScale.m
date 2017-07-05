function [Yhat,shift,scale] = BentIdentityScale(Y)
%scale output data to use bent identity function -std<y<std

%get size
[NNo,SSo] = size(Y);

%mean variance
shift = mean(Y);
scale = ones(1,SSo);

%output
Yhat = ( Y - repmat(shift,NNo,1) ) ./ repmat(scale,NNo,1);

end