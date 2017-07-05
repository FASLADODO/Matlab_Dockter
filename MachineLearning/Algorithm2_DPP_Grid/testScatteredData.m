%test rando data

nn = 500;

mu1 = [5,3];
sigma1 = [1,0.1;0.1,1.0];
mu2 = [1,6];
sigma2 = [1,0.1;0.1,1.0];
mu3 = [10,2];
sigma3 = [1,0.1;0.1,1.0];

r1= mvnrnd(mu1,sigma1,nn);
r2= mvnrnd(mu2,sigma2,nn);
r3= mvnrnd(mu3,sigma3,nn);

mu4 = [2,2];
sigma4 = [1,0.1;0.1,1.0];
mu5 = [4,8];
sigma5 = [1,0.1;0.1,1.0];
mu6 = [8,1];
sigma6 = [1,0.1;0.1,1.0];

r4= mvnrnd(mu4,sigma4,nn);
r5= mvnrnd(mu5,sigma5,nn);
r6= mvnrnd(mu6,sigma6,nn);


Data = [r1;r2;r3;r4;r5;r6];
Labels = [ones(nn,1)*-1;ones(nn,1)*-1;ones(nn,1)*-1;ones(nn,1)*1;ones(nn,1)*1;ones(nn,1)*1];

figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')

%% train grid
split = 7;

[Model] = TrainDPPGridWeak(Data,Labels,split);

[NN,SS] = size(Data);

figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
% gscatter3(Data(:,1),Data(:,2),Data(:,3),Labels)
hold on
for ii = 1:size(Model.limits,1)
    minz = reshape(Model.limits(ii,1,:),1,SS);
    maxz = reshape(Model.limits(ii,2,:),1,SS);
    
    meanz = Model.means(ii,:);
    
    %get line order
    line([minz(1,1),minz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([maxz(1,1),maxz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [minz(1,2),minz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [maxz(1,2),maxz(1,2)])
    hold on
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('manual grid data')

%Try to plot surface with density
%plot grid on data
figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
Surface3D(Model.means(:,1),Model.means(:,2),mean(Model.WeightRegion,2));
hold on
for ii = 1:size(Model.limits,1)
    minz = reshape(Model.limits(ii,1,:),1,SS);
    maxz = reshape(Model.limits(ii,2,:),1,SS);
    
    meanz = Model.means(ii,:);
    
    %get line order
    line([minz(1,1),minz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([maxz(1,1),maxz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [minz(1,2),minz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [maxz(1,2),maxz(1,2)])
    hold on
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('DPP Weights')

%% Try it online

[ClassEst,RawStore] = OnlineDPPGridWeak(Data,Model);

figure
gscatter(Data(:,1),Data(:,2),ClassEst(:,1),'rgc')
hold on
for ii = 1:size(Model.limits,1)
    minz = reshape(Model.limits(ii,1,:),1,SS);
    maxz = reshape(Model.limits(ii,2,:),1,SS);
    
    meanz = Model.means(ii,:);
    
    %get line order
    line([minz(1,1),minz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([maxz(1,1),maxz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [minz(1,2),minz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [maxz(1,2),maxz(1,2)])
    hold on
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('Class Est')

%% plot raw sums

figure
scatter(Data(:,1),Data(:,2),10,RawStore)
hold on
for ii = 1:size(Model.limits,1)
    minz = reshape(Model.limits(ii,1,:),1,SS);
    maxz = reshape(Model.limits(ii,2,:),1,SS);
    
    meanz = Model.means(ii,:);
    
    %get line order
    line([minz(1,1),minz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([maxz(1,1),maxz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [minz(1,2),minz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [maxz(1,2),maxz(1,2)])
    hold on
end
hold off
title('Class Weights')
colorbar 
colormap cool


