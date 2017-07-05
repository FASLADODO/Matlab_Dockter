function K = constructKernel(X1,X2,arg,KernelType)
%done using Foundations of Machine Learning - Mohri
%KernelType: 'Polynomial', 'Gaussian', 'Sigmoid'
%arg: (depends on kernel type) dim, sigma, [a,b]

%See also http://www.mathworks.com/matlabcentral/fileexchange/34864-decision-boundary-using-svms

[NN1,SS1] = size(X1);
[NN2,SS2] = size(X2);
if(NN1 ~= NN2 || SS1 ~= SS2)
    error('X1 and X2 must be same size')
end

switch lower(KernelType)
    case {lower('Linear')}
        K = X1' * X2;
    case {lower('Polynomial')}
        dim = arg(1);
        c = arg(2);
        X1 = reshape(X1',numel(X1),1);
        X2 = reshape(X2',numel(X2),1);
        D  = X1' * X2 + c;
        K = D.^dim;
    case {lower('Gaussian')}
        %really just an RBF
        sigma = arg;
        X1 = reshape(X1',numel(X1),1);
        X2 = reshape(X2',numel(X2),1);
        xn = X1-X2;
        D  = xn'*xn;
%         K = exp(-sigma*D  );
        K = exp(-D / (2*sigma^2) );
    case {lower('Sigmoid')}
        a = arg(1);
        b = arg(2);
        X1 = reshape(X1',numel(X1),1);
        X2 = reshape(X2',numel(X2),1);
        D = a* (X1.*X2) + b;
        K = tanh(D);
    otherwise
        error('KernelType does not exist!');
end

