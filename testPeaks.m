%generate data
data = randn(25,1);

%stanky peaks
[pks,locs] = findpeaks(data)

% plot data with peaks
figure
plot(data)
hold on
scatter(locs,pks,'ro')
hold off

%% 3D

%generate data
data = randn(25,2);

dis = NormRowWise(data);

%stanky peaks
[pks,locs] = findpeaks(dis)

% plot data with peaks
figure
scatter(data(:,1),data(:,2))
hold on
scatter(data(locs,1),data(locs,2),'ro')
hold off