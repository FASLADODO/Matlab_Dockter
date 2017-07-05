function pdf = kde3d(X,X_All,option,h)
%Compute KDE for a 3d data set
% X is the data set X_all is the complete data set (if mixture), h is the bandwidth [h1 h2 h3] (optional)
% option = 'graph' or 'full'

%Example Use
% dat = [randn(1000,1), randn(1000,1), randn(1000,1)];
% pdf = kde3d(dat);
% figure
% scatter3(dat(:,1),dat(:,2),dat(:,3),8,pdf.prob);
% colormap(cool);
% colorbar;
% title('Density Estimates')
% xlabel('x')
% ylabel('y')
% zlabel('z')

%See here:
%http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/ebooks/html/spm/spmhtmlnode18.html

%kernel type 
ktype = 'epanechnikov';

if(nargin < 3)
    option = 'graph';
end
if(nargin < 2)
    X_All = X;
end

%get data set info
[n,dim]=size(X);
[n_all,dim_all]=size(X_All);
if(dim > n) %check on format should be (3 columns x n rows)
    X = X';%resize
    X_All = X_All';%resize
    [n,dim]=size(X);
    [n_all,dim_all]=size(X_All);
end
if(dim > 3)
   error('too many dimensions')
end
%get limits
min_x = min(X_All);
max_x = max(X_All);

if(nargin < 4)
    %optimal bandwidth if not set
    sigma = std(X);
    display('Bandwidths are:')
    h = 1.06*sigma*n^(-(1/6))
    %h = ( 4/(dim+2) )^(1/(dim+4)) * n^(-(1/dim+4)) * sigma
end

pdf.bw = h;

%steps for output vector
steps = 0.01;

%uniformly spaced data (for looping through and plotting)
xlin1 = linspace(min_x(1),max_x(1),n_all*steps);
xlin2 = linspace(min_x(2),max_x(2),n_all*steps);
xlin3 = linspace(min_x(3),max_x(3),n_all*steps);
pdf.xlin =[xlin1; xlin2; xlin3];

switch option
    case 'graph'
        %For individual cloud plot
        ind = 1;
        for jj = 1:n
            disp(sprintf('Loop %i of %i \n', jj, n ));
            %Loop through all relevant x,y values
            sum = 0;
            for ii = 1:n
                %loop through all training data
                val1 = ( X(ii,1) - X(jj,1) ) / h(1);
                val2 = ( X(ii,2) - X(jj,2) ) / h(2);
                val3 = ( X(ii,3) - X(jj,3) ) / h(3);
                %product kernels
                %sum = sum + ( kernels(val1,ktype)*kernels(val2,ktype) );
                %radial kernels
                sum = sum + ( kernel_radial( [val1,val2,val3] ,ktype) );
            end
            pdf.pos(ind,:) = X(jj,:);
            pdf.prob(ind,1) = ( 1/(n*prod(h)) ) * sum; %for plotting
            ind = ind + 1;
        end
    case 'full'
        %for look up table
        ind = 1;
        for xx1 = xlin1 %loop through first dimension
            disp(sprintf('Loop %i of %i \n', xx1, length(xlin1) ));
            for xx2 = xlin2 %loop through second dimension
                for xx3 = xlin3 %loop through third dimension
                    %Loop through all relevant x,y values
                    sum = 0;
                    for ii = 1:n
                        %loop through all training data
                        val1 = ( X(ii,1) - xx1 ) / h(1);
                        val2 = ( X(ii,2) - xx2 ) / h(2);
                        val3 = ( X(ii,3) - xx3 ) / h(3);
                        %product kernels
                        %sum = sum + ( kernels(val1,ktype)*kernels(val2,ktype) );
                        %radial kernels
                        sum = sum + ( kernel_radial( [val1,val2,val3] ,ktype) );
                    end
                    pdf.pos(ind,:) = [xx1 xx2 xx3];
                    pdf.prob(ind,1) = ( 1/(n*prod(h)) ) * sum;
                    ind = ind + 1;
                end
            end
        end
    otherwise
        %For individual cloud plot
        ind = 1;
        for jj = 1:n
            disp(sprintf('Loop %i of %i \n', jj, n ));
            %Loop through all relevant x,y values
            sum = 0;
            for ii = 1:n
                %loop through all training data
                val1 = ( X(ii,1) - X(jj,1) ) / h(1);
                val2 = ( X(ii,2) - X(jj,2) ) / h(2);
                val3 = ( X(ii,3) - X(jj,3) ) / h(3);
                %product kernels
                %sum = sum + ( kernels(val1,ktype)*kernels(val2,ktype) );
                %radial kernels
                sum = sum + ( kernel_radial( [val1,val2,val3] ,ktype) );
            end
            pdf.pos(ind,:) = X(jj,:);
            pdf.prob(ind,1) = ( 1/(n*prod(h)) ) * sum; %for plotting
            ind = ind + 1;
        end
end

%Scale all probabilities to 1
pdf.prob = pdf.prob./max(pdf.prob);
end