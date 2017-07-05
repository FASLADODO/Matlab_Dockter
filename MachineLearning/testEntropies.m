mu1 = 1;
sigma1 = 10;

mu2 = 3;
sigma2 = 10;

nn = 100;

dat1 = mvnrnd(mu1,sigma1,nn);
dat2 = mvnrnd(mu2,sigma2,nn);


p1 = gaussianprob(dat1, mu1, sigma1);
p2 = gaussianprob(dat2, mu2, sigma2);


[ H_X, H_Y, H_XY ] = JointEntropy(p1,p2)



figure
scatter(dat1,zeros(1,nn),'r.')
hold on
scatter(dat1,p1,'c.')
hold on
scatter(dat2,zeros(1,nn),'b.')
hold on
plot(dat2,p2,'g.')

hold off

%% Plottin joints

[d1g,d2g] = meshgrid(dat1,dat2);

[pxg,pyg] = meshgrid(p1,p2);

allP = pxg(:).*pyg(:);
alld = [d1g(:),d2g(:)];

allP = allP./max(allP);

figure
scatter(dat1,ones(1,nn)*mean(dat2),'r.')
hold on
scatter(ones(1,nn)*mean(dat1),dat2,'b.')
hold on
handle = Surface3D(alld(:,1),alld(:,2),allP,'mesh');
hold off

%%
mu1 = [1,3];
sigma1 = [0.5,0; 0,1];

mu2 = [2,4];
sigma2 = [0.5,0; 0,1];

nn = 100;

d1 = mvnrnd(mu1,sigma1,nn);
d2 = mvnrnd(mu2,sigma2,nn);

figure
scatter(d1(:,1),d1(:,2),'r.')
hold on
scatter(d2(:,1),d2(:,2),'b.')
hold off

gg = 60;
CPoint = d1(gg,:);

P = gaussianProbMVarray(d1,sigma1,mu1);
P2 = gaussianProbMV(CPoint,sigma1,mu1)

P(gg)

% plot 2d joints
Data12 = [d1;d2];

limz = [min(Data12);max(Data12)]

xlin = linspace(limz(1,1),limz(2,1),50);
ylin = linspace(limz(1,2),limz(2,2),50);

%mesh grid
[X,Y] = meshgrid(xlin,ylin);

AllGrid = [X(:),Y(:)];

P1 = gaussianProbMVarray(AllGrid,cov(d1),mean(d1));
P2 = gaussianProbMVarray(AllGrid,cov(d2),mean(d2));

jointp = P1.*P2;


figure
scatter(d1(:,1),d1(:,2),'r.')
hold on
scatter(d2(:,1),d2(:,2),'b.')
hold on
handle = Surface3D(AllGrid(:,1),AllGrid(:,2),jointp,'mesh');
hold off

