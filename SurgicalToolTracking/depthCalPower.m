function z = depthCalPower(A,X)

% z = A(1).*(X(:,1).^A(2)) + A(3).*X(:,2) + A(4).*X(:,3) + A(5);
z = A(1).*(X(:,1).^A(2));