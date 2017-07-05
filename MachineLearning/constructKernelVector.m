function K = constructKernelVector(X1,X2,arg,KernelType)
%apply kernels in constructKernel to each row of a matrix

[NN1,SS1] = size(X1);
[NN2,SS2] = size(X2);

if(NN1 ~= NN2)
   if(NN2 == 1)
       X2 = repmat(X2,NN1,1);
   else
       error('mismatch dimensions')
   end
end

for ii = 1:NN1
   x1temp = X1(ii,:);
   x2temp = X2(ii,:);
   K(ii,:) = constructKernel(x1temp,x2temp,arg,KernelType);
end

end