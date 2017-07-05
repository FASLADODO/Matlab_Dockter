% 1D integral

% make some data
nn = 100;
a = randn(nn,1);

%make the function
func = @(x) sum(exp(- norm(x-a) ) );

%some bounds
xmin = -5;
xmax = 5;
xtemp = -5:0.1:5;

%compute intergal
intval = integral(func,xmin,xmax,'ArrayValued',true);

%get rbf at each point
for ii = 1:length(xtemp);
    prob(ii) = func(xtemp(ii));
end

%for comparison
sumscale = sum(prob)

%make plot
figure
scatter(a,zeros(nn,1),'ro')
hold on
plot(xtemp,prob,'b')
hold off
str = sprintf('1D data with integral = %d', intval);
title(str)