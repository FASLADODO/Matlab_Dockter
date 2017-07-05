function [Yhat,shift,scale] = SigmoidScale(Y)
%scale output data to use sigmoid function 0<y<1

%get size
[NNo,~] = size(Y);

%mean variance
shift = min(Y);
scale = range(Y);

%output (0<y<1)
Yhat = ( Y - repmat(shift,NNo,1) ) ./ repmat(scale,NNo,1);

end