%Test linear gaussians

%See LinearGaussianTrain() and LinearGaussianOnline()

% linear gaussian for cube functions

noiz = 100;
params = [3;2.5;-0.1;-0.05];
range = [1,25];
nn = 5000;

xo = linspace(range(1),range(2),nn)';
data = [ones(nn,1), xo, xo.^2, xo.^3];
y = data*params + randn(nn,1)*noiz;

figure
scatter(xo,y)
title('initial data')

%Train this fancy model
ModelLSG = LinearGaussianTrain(data,y);

%test online
[PL,Y_Est] = LinearGaussianOnline(data,y,ModelLSG);

figure
scatter(xo,y,'b.')
hold on
plot(xo,Y_Est,'r-')
hold on
handle2 = Surface3D(xo,y,PL,'mesh');
hold off
title('linear gaussian fit')



%% All in one
[PL,meany] = LinearGaussian(data,y,data,y);

figure
scatter(xo,y,'b.')
hold on
plot(xo,meany,'r-')
hold on
handle2 = Surface3D(xo,y,PL);
hold off
