%% Test leap motion data for seperability
% 
% D1 = load('2016-8-19_15-32-42_circles.csv');
% D2 = load('2016-8-19_15-34-47_ymotion.csv');

% D1 = load('2016-8-22_15-6-11_expert.csv');
% D2 = load('2016-8-22_15-6-57_novice.csv');

% D1 = load('2016-8-24_12-58-53_confident.csv');
% D2 = load('2016-8-24_13-5-35_timid.csv');

D1 = load('2016-8-30_12-53-20_confident.csv');
D2 = load('2016-8-30_12-53-52_timid.csv');


Data = [D1;D2];

Labels = [ones(length(D1),1)*1;ones(length(D2),1)*2];
% LabelKey = {'circles';'ymotion'};
LabelKey = {'expert';'novice'};


key.col.id = 1;
key.col.time = 2;
key.col.x = 3;
key.col.y = 4;
key.col.z = 5;
key.col.segment = 6;
key.col.dx = 7;
key.col.dy = 8;
key.col.dz = 9;
key.col.velmag = 10;
key.col.velalpha = 11;
key.col.velbeta = 12;
key.col.ddx = 13;
key.col.ddy = 14;
key.col.ddz = 15;
key.col.accmag = 16;
key.col.accalpha = 17;
key.col.accbeta = 18;

key.strings = {'frameid','time','x','y','z','segment','dx','dy','dz','velmag','velalpha','velbeta','ddx','ddy','ddz','accmag','accalpha','accbeta'};


figure
gscatter3(Data(:,key.col.x),Data(:,key.col.y),Data(:,key.col.z),Labels)
title('initial data')


cslist = unique(Labels);

DataCopy = Data;

for cc = 1:length(cslist)
    %get current class data
    idx = find(Labels == cslist(cc) );
    Dtemp = DataCopy(idx,:);
    
    %start time from zero;
    Data(idx,key.col.time) = Data(idx,key.col.time) - Dtemp(1,key.col.time);
    
    %get time step
    dt = mean(diff(Dtemp(:,key.col.time)))
    
    %filter the motion just a lil bit
    xtemp = smooth(Dtemp(:,key.col.x));
    ytemp = smooth(Dtemp(:,key.col.y));
    ztemp = smooth(Dtemp(:,key.col.z));
    
    %compute derivatives
    dx = Calculate_velocity( xtemp, dt, 'holoborodko'); 
    dy = Calculate_velocity( ytemp, dt, 'holoborodko'); 
    dz = Calculate_velocity( ztemp, dt, 'holoborodko');
    
    ddx = Calculate_velocity( dx, dt, 'holoborodko'); 
    ddy = Calculate_velocity( dy, dt, 'holoborodko'); 
    ddz = Calculate_velocity( dz, dt, 'holoborodko');
    
    %magnitude and angle of velocity
    velmag = NormRowWise([dx,dy,dz]);
    [velalpha,velbeta] = PitchYaw3D([dx,dy,dz]);
    
    %magnitude and angle of acceleration
    accmag = NormRowWise([ddx,ddy,ddz]);
    [accalpha,accbeta] = PitchYaw3D([ddx,ddy,ddz]);
    
    %store that shiz back in the main matrix
    Data(idx,key.col.x) = xtemp;
    Data(idx,key.col.y) = ytemp;
    Data(idx,key.col.z) = ztemp;
    Data(idx,key.col.dx) = dx;
    Data(idx,key.col.dy) = dy;
    Data(idx,key.col.dz) = dz;
    Data(idx,key.col.velmag) = velmag;
    Data(idx,key.col.velalpha) = velalpha;
    Data(idx,key.col.velbeta) = velbeta;
    Data(idx,key.col.ddx) = ddx;
    Data(idx,key.col.ddy) = ddy;
    Data(idx,key.col.ddz) = ddz;
    Data(idx,key.col.accmag) = accmag;
    Data(idx,key.col.accalpha) = accalpha;
    Data(idx,key.col.accbeta) = accbeta;
end


figure
gscatter3(Data(:,key.col.x),Data(:,key.col.y),Data(:,key.col.z),Labels)
title('filtered data')
xlabel('x')
ylabel('y')
zlabel('z')

figure
gscatter3(Data(:,key.col.dx),Data(:,key.col.dy),Data(:,key.col.dz),Labels)
title('vel data')
xlabel('dx')
ylabel('dy')
zlabel('dz')

figure
gscatter3(Data(:,key.col.velalpha),Data(:,key.col.velbeta),Data(:,key.col.velmag),Labels)
title('vel  mag/vec data')
xlabel('velalpha')
ylabel('velbeta')
zlabel('velmag')


figure
gscatter3(Data(:,key.col.ddx),Data(:,key.col.ddy),Data(:,key.col.ddz),Labels)
title('acc data')
xlabel('ddx')
ylabel('ddy')
zlabel('ddz')


%% Find best states for Seperability

% testcolumns = [key.col.dx,key.col.dy,key.col.dz,key.col.velalpha,key.col.velbeta,key.col.velmag,key.col.accmag,key.col.ddx,key.col.ddy,key.col.ddz,key.col.accalpha,key.col.accbeta,key.col.accmag];
testcolumns = [key.col.velalpha,key.col.velbeta,key.col.velmag,key.col.accalpha,key.col.accbeta,key.col.accmag];
nkcombos = nchoosekmulti(testcolumns,2:length(testcolumns));


[Variations,BestSep,BestSepClass,Data_All,Diff_All] = RelativeRBFSeperabilityCheck(Data,Labels,testcolumns,key.strings,2,'ploton');
BestSep.bestcollabels
columnrbf = BestSep.bestcolumns;

%'accmag'    'ddx' wins 



%% plot the best columns


Dplot = Data(:,columnrbf);
figure
if(length(columnrbf)<3)
    gscatter(Dplot(:,1),Dplot(:,2),Labels,'br')
else
    gscatter3(Dplot(:,1),Dplot(:,2),Dplot(:,3),Labels,'br')
    zlabel(BestSep.bestcollabels(3))
end
title('phase portrait')
xlabel(BestSep.bestcollabels(1))
ylabel(BestSep.bestcollabels(2))


%% use best states in random forest

DataUse = Data(:,columnrbf);

rng(1); % For reproducibility
Mdl = TreeBagger(50,DataUse,Labels,'OOBPrediction','On',...
    'Method','classification')

% view(Mdl.Trees{1},'Mode','graph') %view it but its huge

%plot loss function
figure;
oobErrorBaggedEnsemble = oobError(Mdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';

labelstr = predict(Mdl,DataUse);
labelestimate = str2num(cell2mat(labelstr));

corr = labelestimate == Labels;

acc = mean(corr)

figure
plot(labelestimate)


%% test model

bounds = DataBounds(DataUse);
Grid = ndimgrid(bounds,100);
labelg= predict(Mdl,Grid);

figure
gscatter(Grid(:,1),Grid(:,2),labelg)
title('model bounds')




%% Train RBF discrim
DataTest = Data(:,columnrbf);

%Train the RBF model
[Fits] = RelativeRBFTrain(DataTest,Labels,3,'ploton');

%try it online
thresh = Fits{end}.ThresholdSep;
[Class,Probability] = RelativeRBFOnline(DataTest,Fits,0.3);

Correct = Class == Labels;
Correct(Class == 0) = [];

accuracy = mean(Correct)

%% Get segments and try classifying those

segthresh = 0.3;

SegData = [];
colors = {'r*';'b+'};
ClassEst = [];
ClassActual = [];

plotsegs = false;

if(plotsegs)
	figure
end
for cc = 1:length(cslist)
    idx = find(Labels == cslist(cc) );
    Dtemp = Data(idx,:);
    %get the total number of segments for this class
    totalsegs = max(Dtemp(:,key.col.segment)) - 1;
    ClassActual = [ClassActual; ones(totalsegs,1)*cc];
    
    %loop through all segments
    for ss = 1:totalsegs
        %get data for current segment
        ids = find(Dtemp(:,key.col.segment) == ss);
        sdata = Dtemp(ids,:);
        SegData{cc}.SegID{ss} = sdata;
        %get only the relevent columns (the seperable ones)
        dataon = sdata(:,columnrbf);
        [LLR,Class,ProbClass] = RelativeRBFOnlineLLR(dataon,Fits,segthresh);
        
        ClassEst = [ClassEst; Class];
        
        if(plotsegs)
            plotd = SegData{cc}.SegID{ss}(:,[key.col.x,key.col.y,key.col.z]);
            clf
            scatter3(plotd(:,1),plotd(:,2),plotd(:,3),colors{cc});
            axis([-80,20,-400,200,100,500])
            view(135, 45);
            xlabel('x')
            ylabel('y')
            zlabel('z')
            pause(1)
        end
    end
end


corr = ClassEst == ClassActual;
acc = mean(corr)


%% Try classifying normally

%Have to make data that is segment based, classify each segment with log
%likelihood

DR = Data(:,columnrbf);

figure
gscatter(DR(:,1),DR(:,2),Labels,'br')
title('inital data')

d1 = DR(Labels == 1,:);
d2 = DR(Labels == 2,:);

%now compute radial basis
gamma = 2;
z1 = RadialBasisFunction(d1,gamma);
z2 = RadialBasisFunction(d2,gamma);

%plot in 3D
figure
gscatter(DR(:,1),DR(:,2),Labels,'br')
hold on
scatter3(d1(:,1),d1(:,2),z1,'r*');
handle1 = Surface3D(d1(:,1),d1(:,2),z1);
hold on
handle2 = Surface3D(d2(:,1),d2(:,2),z2);
hold off
title('initial data')

meas = [d1,z1;d2,z2];

MdlLinear = fitcdiscr(meas,Labels);
labest = predict(MdlLinear,meas);

corr = labest == Labels;
acc = mean(corr)

figure
gscatter(DR(:,1),DR(:,2),labest,'br')
title('est class')






