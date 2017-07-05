function z = xyoffsetfit(A,X)
% give X columns as x (1) and y (2) pixel positions, A as coefficients

 z = A(1) + A(2)*X(:,1) + A(3)*X(:,2) + A(4)*(X(:,1).^2) + A(5)*X(:,1).*X(:,2) + A(6)*(X(:,2).^2) + A(7)*(X(:,1).^3) + A(8)*(X(:,1).^2).*X(:,2) + A(9)*X(:,1).*(X(:,2).^2) + A(10)*(X(:,2).^3);     