%test continuous distributions

%http://www.mathworks.com/help/stats/fitgmdist.html

%% 2D

mu1 = [1 2];
Sigma1 = [2 0; 0 0.5];
mu2 = [-3 -5];
Sigma2 = [1 0;0 1];
mu3 = [6 -6];
Sigma3 = [1 0;0 0.5];
rng(1); % For reproducibility
X = [mvnrnd(mu1,Sigma1,1000);mvnrnd(mu2,Sigma2,1000);mvnrnd(mu3,Sigma3,1000)];


%get gaussian mixture model
GMModel = fitgmdist(X,3);

figure
y = [zeros(1000,1);ones(1000,1);ones(1000,1)*2];
h = gscatter(X(:,1),X(:,2),y);
hold on
ezmesh(@(x1,x2)pdf(GMModel,[x1 x2]),get(gca,{'XLim','YLim'}))
title('{\bf Scatter Plot and Fitted Gaussian Mixture Contours}')
legend(h,'Model 0','Model1')
hold off

%% 3D

nn = 1000;

MU1 = [2 2 2];
SIGMA1 = [2 0 0; 0 .5 0; 0 0 1];

MU2 = [-3 -3 -3];
SIGMA2 = [1 0 0; 0 3 0; 0 0 2];

MU3 = [0 4 10];
SIGMA3 = [2 0 0; 0 1 0; 0 0 .5];

X1 = mvnrnd(MU1,SIGMA1,nn);
X2 = mvnrnd(MU2,SIGMA2,nn);
X3 = mvnrnd(MU3,SIGMA3,nn);

figure
scatter3(X1(:,1),X1(:,2),X1(:,3),10,'r.')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),10,'b.')
hold on
scatter3(X3(:,1),X3(:,2),X3(:,3),10,'g.')
hold off

title('Training Distributions')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Group1','Group2','Group3')

X_All = [X1;X2;X3];

%get gaussian mixture model
GMModel = fitgmdist(X_All,3);

F1 = pdf(GMModel,X_All);

figure
scatter3(X_All(:,1),X_All(:,2),X_All(:,3),8,F1(:));

colormap(cool);
colorbar;
title('Density Estimates')
xlabel('x')
ylabel('y')
zlabel('z')

