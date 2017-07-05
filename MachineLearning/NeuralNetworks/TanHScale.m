function [Yhat,shift,scale] = TanHScale(Y)
%scale output data to use TanH function 0<y<inf

%get size
[NNo,~] = size(Y);

%mean variance
shift = mean(Y);
scale = std(Y);

%output (0<y<inf)
Yhat = ( Y - repmat(shift,NNo,1) ) ./ repmat(scale,NNo,1);

end