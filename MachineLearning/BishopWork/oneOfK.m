function x = oneOfK(N,k)

%returns N x k matrix of 0's and 1's, each row only has a single 1

x = zeros(N,k);

%Get some samples
idx = randi([1,k],N,1);
idB = sub2ind( size(x), [1:N]',idx);

x(idB) = 1;

end