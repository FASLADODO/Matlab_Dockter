% test various dimensionality reduction techniques

% Test relieff

%create some fake data
mu1 = [1,1,1];
sigma1 = [0.5,1,1];

mu2 = [5,1,3];
sigma2 = [0.5,1,1];

nn = 1000;

X1 = [];
X2 = [];
for ii = 1:length(mu1)
    X1 = [X1, randn(nn,1)*sigma1(ii) + mu1(ii)];
    X2 = [X2, randn(nn,1)*sigma2(ii) + mu2(ii)];
end


%Get working data
X = [X1;X2];
Labels = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)

figure
gscatter3(X(:,1),X(:,2),X(:,3),Labels);

% test reliefF
K = 50;
[RANKED,WEIGHT] = relieff(X,Labels,K)
[rankrod,weightsrod] = relieffBasic(X,Labels,K)


%% load ionosphere data for relieff

load ionosphere

YI = MapValues(Y,unique(Y),[1;2]);

K = 20;
[RANKED,WEIGHT] = relieff(X,Y,K);

[rankrod,weightsrod] = relieffBasic(X,YI,K);

figure
bar(weightsrod(rankrod));
xlabel('Predictor rank');
ylabel('Predictor importance weight');
title('weights rod')

figure
bar(WEIGHT(RANKED));
xlabel('Predictor rank');
ylabel('Predictor importance weight');
title('weights builtin')

% plot the important ones
figure
gscatter3(X(:,rankrod(1)),X(:,rankrod(2)),X(:,rankrod(3)),YI)
title('good states rod')
xlabel('rank 1')
ylabel('rank 2')
zlabel('rank 3')

figure
gscatter3(X(:,RANKED(1)),X(:,RANKED(2)),X(:,RANKED(3)),YI)
title('good states builtin')
xlabel('rank 1')
ylabel('rank 2')
zlabel('rank 3')


%% compare classifying
Labels = YI;

DataRod = X(:,[rankrod(1),rankrod(2),rankrod(3)]);
DataBuiltin = X(:,[RANKED(1),RANKED(2),RANKED(3)]);


MdlLinearRod = fitcdiscr(DataRod,Labels);
MdlLinearBuiltin = fitcdiscr(DataBuiltin,Labels);

EstRod = predict(MdlLinearRod,DataRod);
EstBuiltin = predict(MdlLinearBuiltin,DataBuiltin);

corrRod = EstRod == Labels;
corrBuiltin = EstBuiltin == Labels;

accRod = mean(corrRod)
accBuiltin = mean(corrBuiltin)


%% compare correlation coefficients

colsuse = [31,32,33];

DataRod = X(:,rankrod(colsuse));
DataBuiltin = X(:,[RANKED(1),RANKED(2),RANKED(3)]);

figure
gscatter3(DataRod(:,1),DataRod(:,2),DataRod(:,3),YI)
title('use states')
xlabel('rank 1')
ylabel('rank 2')
zlabel('rank 3')

Data1 = DataRod(Labels == 1,:);
Data2 = DataRod(Labels == 2,:);

np = min(length(Data1),length(Data2));

R = corrcoef(Data1(1:np,:),Data2(1:np,:))


%% Test PCA

%create some fake data
mu1 = [1,1,1];
sigma1 = [2,1.5,1];

mu2 = [1,2,3];
sigma2 = [1,1.5,2];

nn = 1000;

X1 = [];
X2 = [];
for ii = 1:length(mu1)
    X1 = [X1, randn(nn,1)*sigma1(ii) + mu1(ii)];
    X2 = [X2, randn(nn,1)*sigma2(ii) + mu2(ii)];
end


%Get working data
X = [X1;X2];
Labels = [ones(nn,1)*1; ones(nn,1)*2]; %(-1/1)

figure
gscatter3(X(:,1),X(:,2),X(:,3),Labels);

[coeff_1,score_1,latent_1] = pca(X1);
coeff_1

[coeff_2,score_2,latent_2] = pca(X2);
coeff_2

[coeff_cross,score_cross,latent_cross] = pca(X);
coeff_cross


%% Test pearson correlation

%https://www.mathworks.com/help/stats/quality-of-life-in-u-s-cities.html

load cities
categories

%box plots
figure()
boxplot(ratings,'orientation','horizontal','labels',categories)

%get perason correlation
CMine = correlationPearson(ratings,ratings);

%ignore diagonal
CMine = CMine - diag(diag(CMine));

besties = max(CMine);

for ii = 1:length(CMine)
   [bestmatch(ii),bestidx(ii)] = max(CMine(:,ii));
   figure
   scatter(ratings(:,ii),ratings(:,bestidx(ii)),'ro');
   xlabel(categories(ii,:))
   ylabel(categories(bestidx(ii),:))
   str = sprintf('cat: %s, corr: %f',categories(ii,:),bestmatch(ii));
   title(str)
end

bestmatch
bestidx


%% Test discirminant pearson correlation


%create some fake data
mu1 = [1,1,1];
sigma1 = [1,1,1];

mu2 = [2,3,4];
sigma2 = [1,2,3];

nn = 1000;

X1 = [];
X2 = [];
for ii = 1:length(mu1)
    X1 = [X1, randn(nn,1)*sigma1(ii) + mu1(ii)];
    X2 = [X2, randn(nn,1)*sigma2(ii) + mu2(ii)];
end

%Get working data
X = [X1;X2];
Labels = [ones(nn,1)*1; ones(nn,1)*2]; %(-1/1)

%plottin yo
figure
gscatter3(X(:,1),X(:,2),X(:,3),Labels);

[DPWEIGHT,RANK] = DiscriminantPearson(X,Labels)

%% try discrim pearson on fisheriris

load fisheriris

%get numeric labels instead of char
Labels = MapValues(species,unique(species),[1:length(unique(species))]);

%Get working data
X = meas;

%plottin yo
figure
gscatter3(X(:,1),X(:,2),X(:,3),Labels);

%rods fancy discriminant
[DPWEIGHT,RANK] = DiscriminantPearson(X,Labels)

%gives 4,3,2 as most discriminant

%plottin da besties
figure
gscatter3(X(:,RANK(1)),X(:,RANK(2)),X(:,RANK(3)),Labels);
DPWEIGHT(RANK)

%compare to relieff
K = 10;
[RANKED,WEIGHT] = relieff(X,Labels,K)


%% Test Canonical Correlation

%create some fake data
mu1 = [1,1,1];
sigma1 = [2,1.5,1];

mu2 = [3,4,3];
sigma2 = [1,1.5,2];

nn = 1000;

X1 = [];
X2 = [];
for ii = 1:length(mu1)
    X1 = [X1, randn(nn,1)*sigma1(ii) + mu1(ii)];
    X2 = [X2, randn(nn,1)*sigma2(ii) + mu2(ii)];
end


%Get working data
X = [X1;X2];
Labels = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)

figure
gscatter3(X(:,1),X(:,2),X(:,3),Labels);
title('original data')

[A,B,r,U,V] = canoncorr(X,Labels);
A
B
r

figure
scatter(U(:,1),V(:,1),'r.');
title('U vs V canoncorr')


figure
scatter(1:length(X),X*A,'r.');
title('mapping canoncorr')

params = pinv(X)*Labels

figure
scatter(1:length(X),X*params,'r.');
title('mapping LS')

