function XF = MovingAverage(X,Size)
%Size = 5 for a moving average of 5
%filters each column of X seperately

%http://www.mathworks.com/help/matlab/ref/filter.html
b = (1/Size)*ones(1,Size);
a = 1;
XF = filter(b,a,X);




end