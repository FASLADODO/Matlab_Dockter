nn = 100;

mu1 = [ 1,2];
mu2 = [ 5,6];
sigma1 = [1, 0.1; 0.1 ,2.1];
sigma2 = [1, 0.1; 0.1 ,1.5];

data1 = mvnrnd(mu1,sigma1,nn);
data2 = mvnrnd(mu2,sigma2,nn);

X = [data1;data2]; %2D data
Y = [ones(nn,1)*0;ones(nn,1) ]; %0's and 1's


%get gaussian mixture model
GMModel1 = fitgmdist(data1,1);
F1 = pdf(GMModel1,data1);

GMModel2 = fitgmdist(data2,1);
F2 = pdf(GMModel2,data2);


[Model] = GaussianClassifierTrain(X,Y);
PCX = GaussianClassifierOnline(X,Model);

figure 
plot(PCX)

figure
gscatter(X(:,1),X(:,2),Y);
hold on
handle = Surface3D(data1(:,1),data1(:,2),F1);
hold off