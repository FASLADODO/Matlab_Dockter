%Testing out various ideas for the DPP grid algorithm

%Create some random data
nn = 200;

m1 =1;
m2 =4;
b1 = 130;
b2 = -1;
noise = 16;

tt = linspace(1,100,nn)';
x1 = tt + randn(nn,1)*noise;
x2 = tt + randn(nn,1)*noise;
y1 = x1*m1 + b1 + randn(nn,1)*noise;
y2 = x2*m2 + b2 + randn(nn,1)*noise;

Data = [x1,y1;x2,y2];
% Data = [x1,y1,tt;x2,y2,tt];
Labels = [ones(nn,1)*1;ones(nn,1)*2];

figure
gscatter(Data(:,1),Data(:,2),Labels)
% gscatter3(Data(:,1),Data(:,2),Data(:,3),Labels)
title('sample data')


split = 10;
[NN,SS] = size(Data);

%weaponizes the training
[Model] = TrainDPPGrid(Data,Labels,split);

%plot grid on data
figure
gscatter(Data(:,1),Data(:,2),Labels)
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
gscatter(Data(:,1),Data(:,2),Labels)
hold on
Surface3D(Model.means(:,1),Model.means(:,2),Model.DensityRegion);
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
title('DPP Density')

%Try to plot surface with probs
%plot grid on data
figure
gscatter(Data(:,1),Data(:,2),Labels)
hold on
Surface3D(Model.means(:,1),Model.means(:,2),Model.KLRegion(:,1));
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
title('DPP KL')






