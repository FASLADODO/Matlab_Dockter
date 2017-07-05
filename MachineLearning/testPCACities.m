% test PCA cities data

%load it up
load cities
categories

figure()
boxplot(ratings,'orientation','horizontal','labels',categories)

%% Test the correlation

C = corr(ratings,ratings)

%% get the PCA values built in

w = 1./var(ratings);
[wcoeff,score,latent,tsquared,explained] = pca(ratings,...
'Algorithm','eig','VariableWeights',w);

wcoeff(:,1:3)

%% compare with rods version

[coeff,score,latent,X_PCA] = PCAbasic(ratings,w);

coeff(:,1:3)

%% get orthornormal ratings

coefforth = inv(diag(std(ratings)))*coeff;

%check orthornormal
I = coefforth'*coefforth;
I(1:3,1:3)

% plot scores
figure()
plot(score(:,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

%% classify based on crime

economics = ratings(:,4);
mueconomics = mean(economics);

labels = economics < mueconomics;

Y_pca = ratings*coefforth';

figure()
gscatter(score(:,1),score(:,2),labels)
