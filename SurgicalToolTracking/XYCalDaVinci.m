function z = XYCalDaVinci(A,X)
% give Xcolumns as X (or Y) pixel position (from center) and Z depth

z = A(1)*X(:,1).*X(:,2) + A(2);