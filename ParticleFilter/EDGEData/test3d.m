%create data matrix
nn = 1000;
d = 3;
mu1 = [3 4 3];
scale1 = [2 2 1];
mu2 = [8 6 3];
scale2 = [2 3 2];
mu3 = [-2 8 -1];
scale3 = [2 1 2];

%create data and plot
dat1 = [randn(nn,1)*(scale1(1)/2) + mu1(1), randn(nn,1)*(scale1(2)/2) + mu1(2), randn(nn,1)*(scale1(3)/2) + mu1(3)];
dat2 = [randn(nn,1)*(scale2(1)/2) + mu2(1), randn(nn,1)*(scale2(2)/2) + mu2(2), randn(nn,1)*(scale2(3)/2) + mu2(3)];
dat3 = [randn(nn,1)*(scale3(1)/2) + mu3(1), randn(nn,1)*(scale3(2)/2) + mu3(2), randn(nn,1)*(scale3(3)/2) + mu3(3)];
dat = [dat1;dat2;dat3];
figure
scatter3(dat1(:,1),dat1(:,2),dat1(:,3),'b.')
hold on
scatter3(dat2(:,1),dat2(:,2),dat2(:,3),'r.')
hold on
scatter3(dat3(:,1),dat3(:,2),dat3(:,3),'g.')
hold off
title('original data set')
xlabel('pos')
ylabel('velocity')
zlabel('acceleration')
legend('group1','group2','group3')

% construct the input for BW estimation
pdfh.Mu = dat1';
N = length(dat1);
pdfh.Cov{N} = {} ;
for i = 1 : N
    pdfh.Cov{i} = zeros(d,d) ;
end
pdfh.w = ones(1,N)/N ;
pdfh.n = N;

% estimate the bandwidth
H = estimateBandwidth( pdfh ) ;

disp('The estimated bandwidth:')
H


% compute distro
bw = 0.8;
pdf1 = kde3d(dat1,dat,'graph');
pdf2 = kde3d(dat2,dat,'graph');
pdf3 = kde3d(dat3,dat,'graph');


%%
figure
scatter3(dat1(:,1),dat1(:,2),dat1(:,3),8,pdf1.prob);
hold on
scatter3(dat2(:,1),dat2(:,2),dat2(:,3),8,pdf2.prob);
hold on
scatter3(dat3(:,1),dat3(:,2),dat3(:,3),8,pdf3.prob);
hold off
colormap(cool);
colorbar;
title('Density Estimates')
xlabel('pos')
ylabel('velocity')
zlabel('acceleration')

%%
% compute distro
bw = 0.8;
pdf1 = kde3d(dat1,dat,'full',[0.8,0.8,0.8]);
% pdf2 = kde3d(dat2,dat,'full');
% pdf3 = kde3d(dat3,dat,'full');

figure
scatter3(pdf1.pos(:,1),pdf1.pos(:,2),pdf1.pos(:,3),8,pdf1.prob);
% hold on
% scatter3(dat2(:,1),dat2(:,2),dat2(:,3),8,pdf2.prob);
% hold on
% scatter3(dat3(:,1),dat3(:,2),dat3(:,3),8,pdf3.prob);
% hold off
colormap(cool);
colorbar;
title('Density Estimates')
xlabel('pos')
ylabel('velocity')
zlabel('acceleration')