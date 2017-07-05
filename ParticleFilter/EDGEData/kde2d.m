function pdf = kde2d(X,h)
%X is the data set, h is the bandwidth [h1 h2] (optional)
% Example use
% dat = [randn(1000,1),randn(1000,1)];
% pdf = kde2d(dat);
% [X,Y] = meshgrid(pdf.xlin1,pdf.xlin2);
% f = scatteredInterpolant(pdf.pos(:,1),pdf.pos(:,2),pdf.prob);
% Z = f(X,Y);
% %Plot the interpolated and the nonuniform data to produce:
% figure
% mesh(X,Y,Z) %interpolated
% axis tight; hold on
% scatter3(dat1(:,1),dat1(:,2),zeros(length(dat1),1),'b.')
% hold on
% scatter3(dat2(:,1),dat2(:,2),zeros(length(dat2),1),'r.')
% hold off

%See here:
%http://sfb649.wiwi.hu-berlin.de/fedc_homepage/xplore/ebooks/html/spm/spmhtmlnode18.html

ktype = 'epanechnikov';

%get data set info
[n,dim]=size(X);
if(dim > n) %check on format should be (2 columns x n rows)
    X = X';%resize
    [n,dim]=size(X);
end
if(dim > 2)
   error('too many dimensions')
end
%get limits
min_x = min(X);
max_x = max(X);

if(nargin < 2)
    %optimal bandwidth if not set
    sigma = std(X);
    h = 1.06*sigma*n^(-(1/6))
    %h = ( 4/(dim+2) )^(1/(dim+4)) * n^(-(1/dim+4)) * sigma %bad
end

%steps for output vector
steps = 0.1;

%uniformly spaced data (for looping through and plotting)
xlin1 = linspace(min_x(1),max_x(1),n*steps);
xlin2 = linspace(min_x(2),max_x(2),n*steps);
pdf.xlin1 = xlin1;
pdf.xlin2 = xlin2;

ind = 1;
for xx1 = xlin1 %loop through first dimension
    for xx2 = xlin2 %loop through second dimension
        %Loop through all relevant x,y values
        sum = 0;
        for ii = 1:n
            %loop through all training data
            val1 = ( X(ii,1) - xx1 ) / h(1);
            val2 = ( X(ii,2) - xx2 ) / h(2);
            %product kernels
            %sum = sum + ( kernels(val1,ktype)*kernels(val2,ktype) );
            %radial kernels
            sum = sum + ( kernel_radial( [val1,val2] ,ktype) );
        end
        pdf.pos(ind,:) = [xx1 xx2];
        pdf.prob(ind,1) = ( 1/(n*prod(h)) ) * sum;
        ind = ind + 1;
    end
end

end