function [Yhat,shift,scale] = ReLUScale(Y)
%scale output data to use Relu function 0<y<inf

%get size
[NNo,SSo] = size(Y);

%mean variance
shift = min(Y);
scale = ones(1,SSo);

%output (0<y<inf)
Yhat = ( Y - repmat(shift,NNo,1) ) ./ repmat(scale,NNo,1);

end