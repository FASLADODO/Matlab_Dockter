%% Test PCA

%http://www.cs.otago.ac.nz/cosc453/student_tutorials/principal_components.pdf

%create some fake data
mu1 = [1,2];
sigma1 = [2,1; 1,2];

mu2 = [2,4];
sigma2 = [2,1; 1,2];

nn = 100;

X0 = mvnrnd(mu1,sigma1,nn);
X1 = mvnrnd(mu2,sigma2,nn);


%Get working data
X = [X0;X1];
Y = [ones(nn,1)*-1; ones(nn,1)*1]; %(-1/1)


figure
scatter(X0(:,1),X0(:,2),'r.')
hold on
scatter(X1(:,1),X1(:,2),'b.')
hold off



%% Rods version


%Number of components
NC = 1;

[FeatureVector,X_PCA] = PCAbasic(X,NC);

X0_t = (FeatureVector*X0' )';
X1_t = (FeatureVector*X1' )';

% figure
% scatter(X0_t(:,1),X0_t(:,2),'r.')
% hold on
% scatter(X1_t(:,1),X1_t(:,2),'b.')
% hold off

figure
scatter(X0_t(:,1),ones(nn,1),'r.')
hold on
scatter(X1_t(:,1),ones(nn,1),'b.')
hold off


%% compare with matlab built in

[coeff,score,latent,tsquared] = pca(X);
 
coeff

figure
scatter(X(:,1),X(:,2),20,tsquared)
colormap cool
colorbar







