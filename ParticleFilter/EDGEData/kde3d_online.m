function prob = kde3d_online(X,online)
%Compute KDE for a single point uding 3d training data set
% X is the training data set 
% online is a single point complete data set (if mixture)
% subsample = percentage of data points to use: eg 0.75

%Example Use
% dat = [randn(1000,1), randn(1000,1), randn(1000,1)];
% on = randn(1,3)
% pdf = kde3d(dat,on,0.75);


%See here:
%http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/ebooks/html/spm/spmhtmlnode18.html

%kernel type 
ktype = 'epanechnikov';

%get data set info
[n,dim]=size(X);
if(dim > n) %check on format should be (3 columns x n rows)
    X = X';%resize
    [n,dim]=size(X);
end
if(dim > 3)
   error('too many dimensions')
end


%optimal bandwidth if not set
sigma = std(X);
%display('Bandwidths are:')
h = 1.06*sigma*n^(-(1/6));

sum = 0;
for ii = 1:n
    %loop through all training data
    val1 = ( X(ii,1) - online(1) ) / h(1);
    val2 = ( X(ii,2) - online(2) ) / h(2);
    val3 = ( X(ii,3) - online(3) ) / h(3);
    %product kernels
    %sum = sum + ( kernels(val1,ktype)*kernels(val2,ktype) );
    %radial kernels
    sum = sum + ( kernel_radial( [val1,val2,val3] ,ktype) );
end
prob = ( 1/(n*prod(h)) ) * sum; %for plotting

end