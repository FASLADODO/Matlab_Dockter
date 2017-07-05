%create data matrix
mu1 = 3;
scale1 = 2;
mu2 = 7;
scale2 = 2;

nn = 500;
dat = [randn(nn,1)*(scale1/2) + mu1; randn(nn,1)*(scale2/2) + mu2];

% compute distro
bw = 0.8;
pdf = kde1d(dat);

%plot
figure
scatter(dat,zeros(length(dat),1),'r.');
hold on
plot(pdf.pos,pdf.prob,'b')
hold off