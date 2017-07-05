function z = depthCal(A,X)

z = A(1).*exp(A(2).*X(:,1)) + A(3).*X(:,2) + A(4).*X(:,3) + A(5);
