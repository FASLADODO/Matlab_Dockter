% testDPPGrid_Thesis

%plot the various parameters for training data for the thesis methods

fsize = 14;

nn = 5000;
noisex = 0.1;
noisey = 0.2;

x1 = linspace(0,10,nn)' + randn(nn,1)*noisex;
x2 = linspace(0,10,nn)' + randn(nn,1)*noisex;
d1 = [x1,x1.^2,x1.^0]; %x,x^2,1
d2 = [x2,x2.^2,x2.^0]; %x,x^2,1

%params
p1 = [2.1;-0.09;4];
p2 = [2;-0.12;4];

%output
y1 = d1*p1  + randn(nn,1)*noisey;
y2 = d2*p2  + randn(nn,1)*noisey;

%Make data and labels
Data = [x1,y1; x1,y2];
Labels = [ones(nn,1)*-1; ones(nn,1)*1];

figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthWest';

%% Now make the grid

colormapnew = flipud(cool);

split = 11;

[NN,SS] = size(Data);

[Model] = TrainDPPGrid(Data,Labels,split);

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
    %scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthWest';
% title('manual grid data')

%Try to plot KL
figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
limits = [min(Data); max(Data)]'
Surface3D(Model.means(:,1),Model.means(:,2),Model.KLRegion,'mesh',limits);
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
    %scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
zlabel('W_{KL}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';
% title('DPP KL Divergence')
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'W_{KL}','FontSize',fsize)



%compute S_{j,k} for kicks
% S = (Model.KLRegion .* Model.DensityRegion) ./ (abs(Model.KLRegion) + abs(Model.DensityRegion)) ;
figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
limits = [min(Data); max(Data)]'
Surface3D(Model.means(:,1),Model.means(:,2),Model.SWeightRegion,'mesh',limits);
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
    %scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
zlabel('S_{j,k}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthWest';
%title('DPP KL Divergence')
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'S_{j,k}','FontSize',fsize)

%% classify some nonsense online

nnt = 100;
noiseon = 0.0001;

xon = linspace(0,10,nnt)' + randn(nnt,1)*noiseon;
don = [xon,xon.^2,xon.^0]; %x,x^2,1

%output
yon = don*p1 + randn(nnt,1)*noiseon;

%Make data and labels
Data_Test = [xon,yon];
Labels_Test = [ones(nnt,1)*-1];

%Classify it 
[Class,Score,ScoreTime] = OnlineDPPGrid(Data_Test,Model);

figure
scatter(1:length(ScoreTime),movmean(ScoreTime,3),'r')
xlabel('t','FontSize',fsize)
ylabel('M_{online}','FontSize',fsize)


figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
scatter(Data_Test(:,1),Data_Test(:,2),'k.')
hold off
xlabel('x1','FontSize',fsize)
ylabel('x2','FontSize',fsize)
lh=legend('-1','1','D_{on}');
lh.FontSize = 16;
lh.Location = 'NorthWest';
