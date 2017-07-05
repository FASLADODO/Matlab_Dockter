
nn = 100;

mu = [ 2,3];
sigma = [2, 0; 0 ,3];

data = mvnrnd(mu,sigma,nn);


sigmat = cov(data);
mut = mean(data);



PT = gaussianScaleArray(data,sigmat,mut);


figure
scatter(data(:,1),data(:,2),'r.')
hold on
Surface3D(data(:,1),data(:,2),PT);
hold off