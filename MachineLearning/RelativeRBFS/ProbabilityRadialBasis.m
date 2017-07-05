%% Test Tims RBFS

%% Make some 1D rand data

means{1} = [2];
sigmas{1} = [1];

means{2} = [6];
sigmas{2} = [2];

nn = [1000,1000];

[Data,Labels] = CreateData(means,sigmas,nn);

[Fits] = RelativeRBFTrain(Data,Labels,2,'ploton');

%%

D1 = Data(Labels == 1,:);
D2 = Data(Labels == 2,:);

figure
gscatter(Data,zeros(length(Data),1),Labels,'br');
title('initial data')

%optimal bandwidth (rule of thumb)
sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))
% bw = 2;

cslist = unique(Labels);

for cc = 1:length(cslist)
    %Get class specific data
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    Data_All{cc} = DON;
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    %probabilities
    Probability{cc} = pwithin;
    Opposite{cc} = pbetween;
%     Difference{cc} = pwithin.^1.*log10(pwithin./(pbetween));
    %KL Divergence
    Difference{cc} = (pwithin+pbetween).*log10(pwithin./(pbetween));
end

figure
gscatter(Data,zeros(length(Data),1),Labels,'br');
hold on
scatter(D1,Probability{1},'g.');
hold on
scatter(D2,Probability{2},'c.');
hold off
title('Probability')
xlabel('sample#')
ylabel('P')
ylim([0,1])

figure
gscatter(Data,zeros(length(Data),1),Labels,'br');
hold on
scatter(D1,Difference{1},'g.');
hold on
scatter(D2,Difference{2},'c.');
hold off
title('Seperability')
xlabel('sample#')
ylabel('Relative P')



%% Try classifying with probs


thresh = 2;
mr = sum(cslist); % should be 3 is cslist = [1,2]
%loop through all classes
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    nn = size(D,1);
    %compute both relative probabilities
    DON = Probability{cc};
    DOFF = Opposite{cc};
    COMB = [DON,DOFF];
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    CLASS{cc} = ClassG;
    %record accuracy
    corr = CLASS{cc} == ones(nn,1)*cc;
    AccuracyProb(cc) = mean(corr);
end

AccuracyProb
CLASSALL = [CLASS{1};CLASS{2}];

figure
gscatter(Data,zeros(length(Data),1),CLASSALL,'br');
hold on
scatter(D1,Probability{1},'g.');
hold on
scatter(D2,Probability{2},'c.');
hold off
str = sprintf('Classify w/ Probabilities (%f)',mean(AccuracyProb) );
title(str)
xlabel('sample#')
ylabel('P')

%% Try classifying with discrim


thresh = log10(2);
mr = max(cslist) + 1; % should be 3 is cslist = [1,2]
%loop through all classes
DPlot = [];
ClassPlot = [];
DiffAll = [];
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    
    nn = size(D,1);
    %compute both relative probabilities
    pwithin = Probability{cc};
    pbetween = Opposite{cc};
    
    [PROnline] = ComputeDiscriminantRBFOnline(D,Data_All,bw);

    COMB = [PROnline{cc},PROnline{mr-cc}];
    
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    isgood = (DFIN>thresh);
    CLASS{cc} = isgood.*ClassG;
    %record accuracy
    corr = ClassG == ones(nn,1)*cc;
    AccuracyLimited(cc) = mean(corr(isgood));
    AccuracyAll(cc) = mean(corr);
    
    %extra stuff
    DPlot = [DPlot; D];
    ClassPlot = [ClassPlot;  CLASS{cc}];
    DiffAll = [ DiffAll; PROnline{cc}];
end

AccuracyLimited
AccuracyAll

figure
gscatter(Data,zeros(length(Data),1),ClassPlot,'gbr');
hold on
scatter(D1,Difference{1},'g.');
hold on
scatter(D2,Difference{2},'c.');
hold off
str = sprintf('Classify with Discriminant (%f)',mean(AccuracyLimited) );
title(str)
xlabel('sample#')
ylabel('Relative P')


%% fit to new dist

dcp = linspace(min(DPlot),max(DPlot),500)';

%Scale this
pscale = 1 / max(abs(DiffAll));
DiffAll = DiffAll .*pscale;

good_data = DPlot(DiffAll > thresh,:);
good_data_class = ClassPlot(DiffAll > thresh,:);
good_data_DiffAll = DiffAll(DiffAll > thresh,:);

figure
gscatter(good_data,zeros(length(good_data),1),good_data_class,'br');
hold on
scatter(good_data,good_data_DiffAll,'g.')
hold off

groupings = [1,1];
kmeanlimit = 2;

figure
hold on

Fits = [];
for cc = 1:length(cslist)
    Dfit = good_data(good_data_class == cslist(cc),1);
    
    %T = clusterdata(Dfit,'maxclust',3);
    eva = evalclusters(Dfit,'kmeans','gap','KList',[1:kmeanlimit]);
    groupings(cc) = eva.OptimalK; %assume this is solved
    T = NormalizedGraphCuts(Dfit,groupings(cc));
    idc = 1;
    for ii = 1:max(T)
        Dclust = Dfit(T==ii,:);
        if(length(Dclust) > 10)
            pd = fitdist(Dclust,'Normal');
            Fits{cc}.cluster{idc} = pd;
            idc = idc + 1;
            
            %for plots
            %dcp = sort(Dclust);
            yclust = pdf(pd,dcp);
            plot(dcp,yclust,'LineWidth',2)
            hold on
        end
    end
end
hold on
gscatter(good_data,zeros(length(good_data),1),good_data_class,'br');
hold off
title('new distributions')
xlabel('x data')
ylabel('prob')


%% Make multi modal 1D data

nn = [1000,1000];

means{1} = [2];
sigmas{1} = [1];
means{2} = [6];
sigmas{2} = [1];

[Data1,Labels1] = CreateData(means,sigmas,nn);


means{1} = [8];
sigmas{1} = [1];
means{2} = [-4];
sigmas{2} = [1];

[Data2,Labels2] = CreateData(means,sigmas,nn);

Data = [Data1;Data2];
Labels = [Labels1;Labels2];

plotheight = [ones(nn(1),1)*1; ones(nn(1),1)*1.2;ones(nn(1),1)*1;ones(nn(1),1)*1.2];

D1 = Data(Labels == 1,:);
D2 = Data(Labels == 2,:);

figure
gscatter(Data,plotheight,Labels,'br');
title('initial data')
ylim([0,10])

%optimal bandwidth (rule of thumb)
sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))

cslist = unique(Labels);

for cc = 1:length(cslist)
    %Get class specific data
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    Data_All{cc} = DON;
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    %probabilities
    Probability{cc} = pwithin;
    Opposite{cc} = pbetween;
%     Difference{cc} = pwithin.^1.*log10(pwithin./(pbetween));
    %KL Divergence
    Difference{cc} = (pwithin+pbetween).*log10(pwithin./(pbetween));
end

figure
gscatter(Data,zeros(length(Data),1),Labels,'br');
hold on
scatter(D1,Probability{1},'g.');
hold on
scatter(D2,Probability{2},'c.');
hold off
title('Probability')
xlabel('sample#')
ylabel('P')


figure
gscatter(Data,zeros(length(Data),1),Labels,'br');
hold on
scatter(D1,Difference{1},'g.');
hold on
scatter(D2,Difference{2},'c.');
hold off
title('Seperability')
xlabel('sample#')
ylabel('Relative P')



%% Try classifying with probs


thresh = 2;
mr = sum(cslist); % should be 3 is cslist = [1,2]
%loop through all classes
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    nn = size(D,1);
    %compute both relative probabilities
    DON = Probability{cc};
    DOFF = Opposite{cc};
    COMB = [DON,DOFF];
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    CLASS{cc} = ClassG;
    %record accuracy
    corr = CLASS{cc} == ones(nn,1)*cc;
    AccuracyProb(cc) = mean(corr);
end

AccuracyProb
CLASSALL = [CLASS{1};CLASS{2}];

figure
gscatter(Data,zeros(length(Data),1),CLASSALL,'br');
hold on
scatter(D1,Probability{1},'g.');
hold on
scatter(D2,Probability{2},'c.');
hold off
str = sprintf('Classify w/ Probabilities (%f)',mean(AccuracyProb) );
title(str)
xlabel('sample#')
ylabel('P')

%% Try classifying with discrim


thresh = log10(2);
mr = max(cslist) + 1; % should be 3 is cslist = [1,2]
%loop through all classes
DPlot = [];
ClassPlot = [];
DiffAll = [];
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    
    nn = size(D,1);
    %compute both relative probabilities
    pwithin = Probability{cc};
    pbetween = Opposite{cc};
    
    [PROnline] = ComputeDiscriminantRBFOnline(D,Data_All,bw);

    COMB = [PROnline{cc},PROnline{mr-cc}];
    
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    isgood = (DFIN>thresh);
    CLASS{cc} = isgood.*ClassG;
    %record accuracy
    corr = ClassG == ones(nn,1)*cc;
    AccuracyLimited(cc) = mean(corr(isgood));
    AccuracyAll(cc) = mean(corr);
    
    %extra stuff
    DPlot = [DPlot; D];
    ClassPlot = [ClassPlot;  CLASS{cc}];
    DiffAll = [ DiffAll; PROnline{cc}];
end

AccuracyLimited
AccuracyAll

figure
gscatter(Data,zeros(length(Data),1),ClassPlot,'gbr');
hold on
scatter(D1,Difference{1},'g.');
hold on
scatter(D2,Difference{2},'c.');
hold off
str = sprintf('Classify with Discriminant (%f)',mean(AccuracyLimited) );
title(str)
xlabel('sample#')
ylabel('Relative P')


%% fit to new dist

dcp = linspace(min(DPlot),max(DPlot),500)';

good_data = DPlot(DiffAll > thresh,:);
good_data_class = ClassPlot(DiffAll > thresh,:);

figure
gscatter(good_data,zeros(length(good_data),1),good_data_class,'br');

groupings = [2,2];

figure
hold on

Fits = [];
for cc = 1:length(cslist)
    Dfit = good_data(good_data_class == cslist(cc),1);
    
    %T = clusterdata(Dfit,'maxclust',3);
    T = NormalizedGraphCuts(Dfit,groupings(cc));
    idc = 1;
    for ii = 1:max(T)
        Dclust = Dfit(T==ii,:);
        if(length(Dclust) > 10)
            pd = fitdist(Dclust,'Normal');
            Fits{cc}.cluster{idc} = pd;
            idc = idc + 1;
            
            %for plots
            %dcp = sort(Dclust);
            yclust = pdf(pd,dcp);
            plot(dcp,yclust,'LineWidth',2)
            hold on
        end
    end
end
hold on
gscatter(good_data,zeros(length(good_data),1),good_data_class,'br');
hold off
title('new distributions')
xlabel('x data')
ylabel('prob')



%% Now try with 2D

means{1} = [2,2];
sigmas{1} = [1,0;0,1];

means{2} = [4,4];
sigmas{2} = [2,0;0,2];

nn = [500,500];

[Data,Labels] = CreateData(means,sigmas,nn);

D1 = Data(Labels == 1,:);
D2 = Data(Labels == 2,:);

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
title('initial data')

%optimal bandwidth (rule of thumb)
sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))
% bw = 2;

cslist = unique(Labels);

for cc = 1:length(cslist)
    %Get class specific data
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    DataC{cc} = DON;
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    %probabilities
    Probability{cc} = pwithin;
    Opposite{cc} = pbetween;
%     Difference{cc} = pwithin.^1.*log10(pwithin./(pbetween));
    %KL Divergence
    Difference{cc} = (pwithin+pbetween).*log10(pwithin./(pbetween));
end

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
cs = 1;
h1 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
hold on
cs = 2;
h2 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
hold off
title('Probability')
xlabel('sample x')
ylabel('sample y')
zlabel('P')
zlim([0,1])
view(45, 20);

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
cs = 1;
h1 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Opposite{cs});
hold on
cs = 2;
h2 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Opposite{cs});
hold off
title('Opposite')
xlabel('sample x')
ylabel('sample y')
zlabel('P')
view(45, 20);

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
cs = 1;
h1 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Difference{cs});
hold on
cs = 2;
h2 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Difference{cs});
hold off
title('Seperability')
xlabel('sample x')
ylabel('sample y')
zlabel('Relative P')
%zlim([0,100])
view(45, 20);

%% Try classifying with probs


thresh = 2;
mr = sum(cslist); % should be 3 is cslist = [1,2]
%loop through all classes
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    nn = size(D,1);
    %compute both relative probabilities
    DON = Probability{cc};
    DOFF = Opposite{cc};
    COMB = [DON,DOFF];
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    CLASS{cc} = ClassG;
    %record accuracy
    corr = CLASS{cc} == ones(nn,1)*cc;
    AccuracyProb(cc) = mean(corr);
end

AccuracyProb
CLASSALL = [CLASS{1};CLASS{2}];

figure
gscatter(Data(:,1),Data(:,2),CLASSALL,'br');
hold on
cs = 1;
h1 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
hold on
cs = 2;
h2 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
hold off
str = sprintf('Classify w/ Probabilities (%f)',mean(AccuracyProb) );
title(str)
xlabel('sample x')
ylabel('sample y')
zlabel('P')
view(45, 20);

%% Try classifying with discrim


thresh = 0.1;
mr = max(cslist) + 1; % should be 3 is cslist = [1,2]
%loop through all classes
DPlot = [];
SPlot = [];
ClassPlot = [];
DiffAll = [];
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    nn = size(D,1);
    %compute both relative probabilities
    pwithin = Probability{cc};
    pbetween = Opposite{cc};
    
    [PROnline] = ComputeDiscriminantRBFOnline(D,DataC,bw);
    
    COMB = [PROnline{cc},PROnline{mr-cc}];
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    isgood = (DFIN>thresh);
    CLASS{cc} = isgood.*ClassG;
    %record accuracy
    corr = ClassG == ones(nn,1)*cc;
    AccuracyLimited(cc) = mean(corr(isgood));
    AccuracyAll(cc) = mean(corr);
    
    %extra stuff
    DPlot = [DPlot; D];
    ClassPlot = [ClassPlot;  CLASS{cc}];
    DiffAll = [ DiffAll; PROnline{cc}];
    SPlot{cc} = PROnline{cc};
end

AccuracyLimited
AccuracyAll
SEPALL = [SPlot{1};SPlot{2}];


figure
gscatter(Data(:,1),Data(:,2),ClassPlot,'gbr');
hold on
cs = 1;
h1 = Surface3D(D1(:,1),D1(:,2),SPlot{cs});
hold on
cs = 2;
h2 = Surface3D(D2(:,1),D2(:,2),SPlot{cs});
hold off
str = sprintf('Classify w/ Discriminant (%f)',mean(AccuracyLimited) );
title(str)
xlabel('sample x')
ylabel('sample y')
zlabel('Relative P')
% zlim([0,100])
view(45, 20);


%% fit to new dist 2-D


good_data = DPlot(DiffAll > thresh,:);
good_data_class = ClassPlot(DiffAll > thresh,:);

bounds = DataBounds(good_data);
dcp = ndimgrid(bounds,500);

groupings = [1,1];

figure

Fits = [];
for cc = 1:length(cslist)
    Dfit = good_data(good_data_class == cslist(cc),:);
    
    %T = clusterdata(Dfit,'maxclust',3);
    T = NormalizedGraphCuts(Dfit,groupings(cc));
    idc = 1;
    for ii = 1:max(T)
        Dclust = Dfit(T==ii,:);
        if(length(Dclust) > 10)
            mu = mean(Dclust);
            sigma = cov(Dclust);
            Fits{cc}.cluster{idc}.mean = mu;
            Fits{cc}.cluster{idc}.sigma = sigma;
            Fits{cc}.cluster{idc}.points = length(Dclust);
            idc = idc + 1;
            
            %for plots
            %dcp = Dclust;
            %yclust = pdf(pd,dcp);
            yclust = mvnpdf(dcp,mu,sigma);
            Surface3D(dcp(:,1),dcp(:,2),yclust);
            hold on
        end
    end
end
hold on
gscatter(good_data(:,1),good_data(:,2),good_data_class,'br');
hold off
title('new distributions')
xlabel('x data')
ylabel('prob')

%% Phase portrait data

fontS = 14;

paramNoise = 0.1; %0.2, 0.7
statenoise = 0.1;
runs = 100;
nn = 100;
x = linspace(0,2.5,nn)';

paramsinit(1,:) = [-2,1.7];
paramsinit(2,:) = [-2.4,2.0];

Data = [];
Labels = [];
for rr = 1:runs
    xon = x + randn(nn,1)*statenoise;
    
    %parameter noise
    for jj = 1:size(paramsinit,1)
        for ii = 1:length(paramsinit(jj,:))
            params(jj,ii) = paramsinit(jj,ii) + (rand(1)-.5)*paramNoise*max(paramsinit(jj,:));
        end
    end

    %Nonlinear equation
    xdot1 = -params(1,1) + params(1,1).*exp(-params(1,2).*xon);
    xdot2 = -params(2,1) + params(2,1).*exp(-params(2,2).*xon);
    
    %state noise
    xdot1 = xdot1 + randn(nn,1)*statenoise;
    xdot2 = xdot2 + randn(nn,1)*statenoise;

    Data = [Data; xon,xdot1; xon,xdot2];
    Labels = [Labels; ones(nn,1)*1; ones(nn,1)*2]; 
end

D1 = Data(Labels == 1,:);
D2 = Data(Labels == 2,:);

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
title('initial data')

%optimal bandwidth (rule of thumb)
sig = norm(std(Data));
n = length(Data);
bw = 1.06*sig*(n^(-1/5))
gamma = 1/(2*bw^2)
%bw = 1;

cslist = unique(Labels);

for cc = 1:length(cslist)
    %Get class specific data
    DON = Data(Labels == cslist(cc),:);
    DOFF = Data(Labels ~= cslist(cc),:);
    DataC{cc} = DON;
    %compute dem rbfs
    [pwithin,pbetween] = ComputeDiscriminateRBF(DON,DOFF,bw);
    %probabilities
    Probability{cc} = pwithin;
    Opposite{cc} = pbetween;
%     Difference{cc} = pwithin.^1.*log10(pwithin./(pbetween));
    %KL Divergence
    Difference{cc} = (pwithin+pbetween).*log10(pwithin./(pbetween));
end

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
% cs = 1;
% h1 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
% hold on
% cs = 2;
% h2 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
h1 = scatter3(D1(:,1),D1(:,2),Opposite{1},'c.');
hold on
h2 = scatter3(D2(:,1),D2(:,2),Opposite{2},'g.');
hold off
title('Opposite')
xlabel('sample x')
ylabel('sample y')
zlabel('P')
view(45, 20);


figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
h1 = scatter3(D1(:,1),D1(:,2),Probability{1},'c.');
hold on
h2 = scatter3(D2(:,1),D2(:,2),Probability{2},'g.');
hold off
title('Probability')
xlabel('sample x')
ylabel('sample y')
zlabel('P')
view(45, 20);

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
% cs = 1;
% h1 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Difference{cs});
% hold on
% cs = 2;
% h2 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Difference{cs});
h1 = scatter3(D1(:,1),D1(:,2),Difference{1},'c.');
hold on
h2 = scatter3(D2(:,1),D2(:,2),Difference{2},'g.');
hold off
title('Seperability')
xlabel('sample x')
ylabel('sample y')
zlabel('Relative P')
%zlim([0,100])
view(45, 20);


%% Try classifying with probs

mr = sum(cslist); % should be 3 is cslist = [1,2]
DPlot = [];
%loop through all classes
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    DPlot = [DPlot; D];
    nn = size(D,1);
    %compute both relative probabilities
    PON = Probability{cc};
    POFF = Opposite{cc};
    COMB = [PON,POFF];
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    CLASS{cc} = ClassG;
    %record accuracy
    corr = CLASS{cc} == ones(nn,1)*cc;
    AccuracyProb(cc) = mean(corr);
end

AccuracyProb
CLASSALL = [CLASS{1};CLASS{2}];

figure
gscatter(DPlot(:,1),DPlot(:,2),CLASSALL,'br');
hold on
cs = 1;
h1 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
hold on
cs = 2;
h2 = Surface3D(DataC{cs}(:,1),DataC{cs}(:,2),Probability{cs});
hold off
str = sprintf('Classify w/ Probabilities (%f)',mean(AccuracyProb) );
title(str)
xlabel('sample x')
ylabel('sample y')
zlabel('P')
view(45, 20);

%% Try classifying with discrim


thresh = 0.1; %log10(2);
mr = max(cslist) + 1; % should be 3 is cslist = [1,2]
DPlot = [];
SPlot = [];
%loop through all classes
for cc = 1:length(cslist)
    %get the current class data and size
    lkp = [cc,mr-cc];
    D = Data(Labels == cslist(cc),:);
    DPlot = [DPlot; D];
    nn = size(D,1);
    %compute both relative probabilities
    pwithin = Probability{cc};
    pbetween = Opposite{cc};
    
    [PROnline] = ComputeDiscriminantRBFOnline(D,DataC,bw);
    %for make sweet plotz
    SPlot{cc} = PROnline{cc};
    
    COMB = [PROnline{cc},PROnline{mr-cc}];
    %find the max ratio
    [DFIN,DID] = max(COMB,[],2);
    %convert to actual class
    ClassG = lkp(DID)';
    %only use classification if above a threshold
    isgood = (DFIN>thresh);
    CLASS{cc} = isgood.*ClassG;
    %record accuracy
    corr = ClassG == ones(nn,1)*cc;
    AccuracyLimited(cc) = mean(corr(isgood));
    AccuracyAll(cc) = mean(corr);
end

AccuracyLimited
AccuracyAll
CLASSALL = [CLASS{1};CLASS{2}];
SEPALL = [SPlot{1};SPlot{2}];

figure
gscatter(DPlot(:,1),DPlot(:,2),CLASSALL,'gbr');
hold on
cs = 1;
h1 = Surface3D(D1(:,1),D1(:,2),SPlot{cs});
hold on
cs = 2;
h2 = Surface3D(D2(:,1),D2(:,2),SPlot{cs});
hold off
str = sprintf('Classify w/ Discriminant (%f)',mean(AccuracyLimited) );
title(str)
xlabel('sample x')
ylabel('sample y')
zlabel('Relative P')
% zlim([0,100])
view(45, 20);


%% Try normalized graph cuts on the seperable data

AllD = [DPlot,SEPALL];

GoodD = DPlot(SEPALL > thresh,:);


figure
scatter(DPlot(:,1),DPlot(:,2),'b');
hold on
scatter(GoodD(:,1),GoodD(:,2),'r');
hold off

k = length(cslist);
LABYA = NormalizedGraphCuts(GoodD,k);

figure
gscatter(GoodD(:,1),GoodD(:,2),LABYA,'br');


%% Try fitting KL Divergence XY for each class with least squares


Datarz = [DPlot, DPlot(:,1).*DPlot(:,2), ones(length(DPlot),1)];
Yout = SEPALL;

ParamSep = pinv(Datarz)*Yout


[X,Y] = meshgrid(linspace(0,2.5,50)', linspace(0,3,50)');
plottall =[X(:), Y(:), X(:).*Y(:), ones(length(X(:)),1) ];

zplot = plottall*ParamSep;


figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
hold on
h1 = Surface3D(plottall(:,1),plottall(:,2),zplot);
hold off
title('LS seperability')
xlabel('sample x')
ylabel('sample y')
zlabel('classifiability')
%zlim([0,100])
view(45, 20);


