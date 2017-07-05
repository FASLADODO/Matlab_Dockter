%% try it on fisher for paper sake

load fisheriris.mat

X = meas;

csstr = unique(species)';
Labels = MapValues(species,csstr,[1,2,3]);

key = {'PW','PL','SW','SL'};


id = Labels ~= 1;
X = X(id,:);
Labels = Labels(id,:);

figure
gscatter3(X(:,1),X(:,2),X(:,3),csstr(Labels))
xlabel(key(1))
ylabel(key(2))
zlabel(key(3))
title('Fisher Data ')


%% SCALE ALL DIMENSIONS TO HAVE THE SAME RANGE (See Lin2006 Paper)

XRaw = X;

%Scaling function based on variance and mean
X = NormalizeFeatures(XRaw);

figure
gscatter3(XRaw(:,1), XRaw(:,2), XRaw(:,3),csstr(Labels))
xlabel(key(1))
ylabel(key(2))
zlabel(key(3))
title('UnScaled Fisher Data')

figure
gscatter3(X(:,1), X(:,2), X(:,3),csstr(Labels))
xlabel(key(1))
ylabel(key(2))
zlabel(key(3))
title('Scaled Fisher Data')


%% Run relief f on all states cuz why the hell not

%Best
% beststates = 'SL'    'SW'    'PW'    'PL'
% weights = 0.0229    0.0173    0.0061    0.0003

[RANK,WEIGHT] = relieff(X,Labels,7);


beststates = key(RANK)
WEIGHT(RANK)

plotcols = [RANK(1),RANK(2),RANK(3)];
figure
gscatter3(X(:,plotcols(1)),X(:,plotcols(2)),X(:,plotcols(3)),csstr(Labels));
xlabel(key(plotcols(1)))
ylabel(key(plotcols(2)))
zlabel(key(plotcols(3)))
title('top 3 relieff fisher')

%% run reliefrbf with scaled data

%Best:
% weights = [1.9748    0.3350    0.2946]
% overall =0.8681
% column numbers = [1     3     4]
% labels = ['PW'    'SW'    'SL']


testcolumns = [1,2,3,4];
nkcombos = nchoosekmulti(testcolumns,2:length(testcolumns))


[Variations,BestSep,BestSepClass,Data_All,Diff_All] = RelativeRBFSeperabilityCheck(X,Labels,testcolumns,key,1,'ploton');
BestSep.bestcollabels
BestSep.metric
idn = 2;
BestSepClass{idn}.max
BestSepClass{idn}.metric
plotcols = BestSepClass{idn}.bestcolumns
collabs = BestSepClass{idn}.bestcollabels

figure
gscatter3(X(:,plotcols(1)),X(:,plotcols(2)),X(:,plotcols(3)),csstr(Labels));
xlabel(key(plotcols(1)))
ylabel(key(plotcols(2)))
zlabel(key(plotcols(3)))
title('top 3 Scaled-reliefrbf fisher Scaled')

figure
gscatter3(XRaw(:,plotcols(1)),XRaw(:,plotcols(2)),XRaw(:,plotcols(3)),csstr(Labels));
xlabel(key(plotcols(1)))
ylabel(key(plotcols(2)))
zlabel(key(plotcols(3)))
title('top 3 Scaled-reliefrbf fisher Raw')

%% run reliefrbf with raw data

%Best:
% weights = [3.6834    0.4390    0.3552]
% overall = 1.4925
% column numbers = [2     3     4]
% labels = ['PL'    'SW'    'SL']

testcolumns = [1,2,3,4];
nkcombos = nchoosekmulti(testcolumns,2:length(testcolumns))


[Variations,BestSepRaw,BestSepClassRaw,Data_AllRaw,Diff_AllRaw] = RelativeRBFSeperabilityCheck(XRaw,Labels,testcolumns,key,1,'ploton');
BestSepRaw.bestcollabels
BestSepRaw.metric
idn = 3;
BestSepClassRaw{idn}.max
BestSepClassRaw{idn}.metric
plotcolsR = BestSepClassRaw{idn}.bestcolumns
collabs = BestSepClassRaw{idn}.bestcollabels

figure
gscatter3(X(:,plotcolsR(1)),X(:,plotcolsR(2)),X(:,plotcolsR(3)),csstr(Labels));
xlabel(key(plotcolsR(1)))
ylabel(key(plotcolsR(2)))
zlabel(key(plotcols(3)))
title('top 3 raw-reliefrbf fisher Scaled')

figure
gscatter3(XRaw(:,plotcolsR(1)),XRaw(:,plotcolsR(2)),XRaw(:,plotcolsR(3)),csstr(Labels));
xlabel(key(plotcolsR(1)))
ylabel(key(plotcolsR(2)))
zlabel(key(plotcolsR(3)))
title('top 3 raw-reliefrbf fisher Raw')


%%
fsize = 16;

% cc = 2;
% scatter3(BestSepClass{idn}.DataClass{cc}(:,1),BestSepClass{idn}.DataClass{cc}(:,2),BestSepClass{idn}.Difference{cc});

%Plot the seperation for 2 columns with surface
idn =2;
figure
gscatter(XRaw(:,plotcols(1)),XRaw(:,plotcols(2)),csstr(Labels)','gk','x+',8);
hold on
cc = 1;
Surface3D(BestSepClassRaw{idn}.DataClass{cc}(:,1),BestSepClassRaw{idn}.DataClass{cc}(:,2),BestSepClass{idn}.Difference{cc});
hold on
cc = 2;
Surface3D(BestSepClassRaw{idn}.DataClass{cc}(:,1),BestSepClassRaw{idn}.DataClass{cc}(:,2),BestSepClass{idn}.Difference{cc});
% hold on
% cc = 3;
% Surface3D(BestSepClassRaw{idn}.DataClass{cc}(:,1),BestSepClassRaw{idn}.DataClass{cc}(:,2),BestSepClass{idn}.Difference{cc});
hold off
% title('best seperation Fisher Iris')
xlabel(key(plotcols(1)),'FontSize',fsize)
ylabel(key(plotcols(2)),'FontSize',fsize)
zlabel('$W_{rbf}$','FontSize',fsize,'interpreter','latex')
hc = colorbar;
ylabel(hc, 'W_{rbf}','FontSize',fsize)
lgz= findobj(gcf,'tag','legend'); 
set(lgz,'FontSize',fsize)


%Plot the seperation for 2 columns with surface
idn =2;
figure
gscatter(XRaw(:,plotcols(1)),XRaw(:,plotcols(2)),csstr(Labels)','gk','x+',8);
% title('best seperation Fisher Iris')
xlabel(key(plotcols(1)),'FontSize',fsize)
ylabel(key(plotcols(2)),'FontSize',fsize)
lgz= findobj(gcf,'tag','legend'); 
set(lgz,'FontSize',fsize)

