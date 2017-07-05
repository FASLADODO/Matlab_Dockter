%create data matrix
nn = 1000;
mu1 = [3 4];
scale1 = [2 2];
mu2 = [8 6];
scale2 = [2 3];

%create data and plot
dat1 = [randn(nn,1)*(scale1(1)/2) + mu1(1), randn(nn,1)*(scale1(2)/2) + mu1(2)];
dat2 = [randn(nn,1)*(scale2(1)/2) + mu2(1), randn(nn,1)*(scale2(2)/2) + mu2(2)];
dat = [dat1;dat2];
figure
scatter(dat1(:,1),dat1(:,2),'b.')
hold on
scatter(dat2(:,1),dat2(:,2),'r.')
hold off
title('original data set')
xlabel('pos')
ylabel('velocity')

% compute distro
bw = 0.8;
pdf = kde2d(dat);

[X,Y] = meshgrid(pdf.xlin1,pdf.xlin2);
%The key to this process is to use scatteredInterpolant to interpolate the values of the function at the uniformly spaced points, based on the values of the function at the original data points (which are random in this example). This statement uses the default linear interpolation to generate the new data:
f = scatteredInterpolant(pdf.pos(:,1),pdf.pos(:,2),pdf.prob);
Z = f(X,Y);

%Plot the interpolated and the nonuniform data to produce:
figure
mesh(X,Y,Z) %interpolated
axis tight; hold on
scatter3(dat1(:,1),dat1(:,2),zeros(length(dat1),1),'b.')
hold on
scatter3(dat2(:,1),dat2(:,2),zeros(length(dat2),1),'r.')
hold off

colormap cool