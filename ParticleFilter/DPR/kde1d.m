function f = kde1d(X,h)
%X is the data set, h is the bandwidth
%example use
% bw = 0.8;
% dat = randn(1000,1);
% f = kde1d(dat,bw);
% figure
% scatter(dat,zeros(1000,1),'r.');
% hold on
% plot(f.pos,f.prob,'b')
% hold off

%kernel type
ktype = 'epanechnikov';

%get data set info
n=length(X);

if(nargin == 1)
    %optimal bandwidth if not set
    sigma = std(X)
    h = 1.5*sigma*n^(-(1/6))
end

%steps for x vector
mind = min(X);
maxd = max(X);
steps = (maxd - mind)/(n*h);


%loop through new x values and all training data to compute pdf
ind = 1;
for xx = mind:steps:maxd 
    %Loop through all relevant x values
    sum = 0;
    for ii = 1:n
        %loop through all training data
        val = ( X(ii) - xx ) / h; %compute distance to training data scaled by bandwidth
        sum = sum + kernels(val,ktype);
    end
    f.pos(ind) = xx;
    f.prob(ind) = ( 1/(h*n) ) * sum;
    ind = ind + 1;
end

end

