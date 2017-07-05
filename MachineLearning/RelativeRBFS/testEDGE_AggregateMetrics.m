%% first load in this large data structure, if isnt already

%load dataglb
if( ~exist('DataGlb'))
    fprintf('Loading EdgeDataGlb.mat (huge) ...');
    try
        load('C:\temp\EdgeDataGlb.mat')
    catch
        try
            load('D:\temp\EdgeDataGlb.mat')
        catch
            disp('Could not load local. Loading from M drive...');
            load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EdgeDataGlb.mat')
        end
    end
    try
        load('C:\temp\EDGE_Segments.mat')
    catch
        try
            load('D:\temp\EDGE_Segments.mat')
        catch
            disp('Could not load local. Loading from M drive...');
            load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EDGE_Segments.mat')
    
        end
    end
    fprintf(' DONE.\n');
end



% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tsk = 3; %PegTx:1 Cutting:2 Suturing:3
TaskStr = DataGlb.Tasks;

myGroups = [g.flsNov g.flsInt g.gtExp];
groupStr = {'Novice','Intermediate','Expert'};

% Choose Segment Type
%SegType = SegS.betweenFgActs; % (based only on Forces, excludes position
%SegType = SegS.betweenZSpd; % (once the tool slows to zero speed)
SegType = SegS.betweenGrActs; % (includes forces and grasper position)

%Get column names
contentcolumns = DataGlb.content(2,:);
contentDescriptions = DataGlb.content(3,:);

%FLS scores
flsColumn = DataGlb.lookupCol('FLS-score');
timeColumn =  DataGlb.lookupCol('General:TotalTime-homebase');
toolPathColumn = DataGlb.lookupCol('General:TotalToolpath-BothHands');
EOMColumn = DataGlb.lookupCol('General:EoM-BothHands');
% JerkColumn = DataGlb.lookupCol('General: Smoothness Jerk');
JerkColumn = DataGlb.lookupCol('MeanToolJerk-Both');
CrvColumn = DataGlb.lookupCol('MeanToolPathCrv-Both');
AccColumn = DataGlb.lookupCol('MeanToolAcc-Both');

%see description
DataGlb.content(3,JerkColumn);

%test
alljerk = cell2mat(DataGlb.content(DataGlb.LogIdx{tsk},JerkColumn));

%for access
key.c.PathLength = 1;
key.c.EOM = 2;
key.c.Jerk = 3;
key.c.Crv = 4;
key.c.Acc = 5;
key.c.Time = 6; %not used
key.c.All = {'Path Length','EOM','Jerk','Crv','Acc','Task Time'};


%% Transform each segment to a new coordinate system: Origin is the end of the segment

%for storage
DataAll = [];
LabelsAll = [];

idx = 1;

fprintf('Getting Group Data...');
for gg = 1:length(myGroups) %looping through array of logs
    
    
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)})
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        %Get metrics for trial
        PathLength_Trial = cell2mat(DataGlb.content(i,toolPathColumn));
        EOM_Trial = cell2mat(DataGlb.content(i,EOMColumn));
        Jerk_Trial = cell2mat(DataGlb.content(i,JerkColumn));
        Crv_Trial = cell2mat(DataGlb.content(i,CrvColumn));
        Acc_Trial = cell2mat(DataGlb.content(i,AccColumn));
        Time_Trial = cell2mat(DataGlb.content(i,timeColumn));
        
        %stash it
        DataAll(idx,key.c.PathLength) = PathLength_Trial ; %to scale for time effects
        DataAll(idx,key.c.EOM) = EOM_Trial;
        DataAll(idx,key.c.Jerk) = Jerk_Trial;
        DataAll(idx,key.c.Crv) = Crv_Trial;
        DataAll(idx,key.c.Acc) = Acc_Trial;
        DataAll(idx,key.c.Time) = Time_Trial;
        LabelsAll(idx,:) = gg;
        
        %increment
        idx = idx + 1;
    end
end
fprintf('Done \n');

%for plotting
pc = [key.c.PathLength,key.c.EOM,key.c.Acc];

figure
gscatter3(DataAll(:,pc(1)),DataAll(:,pc(2)),DataAll(:,pc(3)),groupStr(LabelsAll));
xlabel(key.c.All(pc(1)))
ylabel(key.c.All(pc(2)))
zlabel(key.c.All(pc(3)))
title('sample state plot')


pc2 = [key.c.Time,key.c.PathLength];

figure
gscatter(DataAll(:,pc2(1)),DataAll(:,pc2(2)),groupStr(LabelsAll)');
xlabel(key.c.All(pc2(1)))
ylabel(key.c.All(pc2(2)))
title('is metric time scaled?')


%% only keep novice and expert

isgroup = ismember(LabelsAll,[1,3]);

%store in new data matrix
Data = DataAll(isgroup,:);
Labels = LabelsAll(isgroup,:);

noviceCount = sum(ismember(LabelsAll,[1]))
expertCount = sum(ismember(LabelsAll,[3]))

%% Now lets to leave one out 

fc = [key.c.PathLength,key.c.EOM,key.c.Crv,key.c.Jerk];
% fc = [key.c.PathLength];

epochs = length(Labels);

stashClassification = [];
stashData = [];

fprintf('Running Leave One Out...');
for ii = 1:epochs
    %only keep training data
   TrainData = Data(:,fc);
   TrainLabel = Labels;
   TrainData(ii,:) = [];
   TrainLabel(ii,:) = [];
   
   %only test on leave out
   TestData = Data(ii,fc);
   TestLabel = Labels(ii,:);
   
   %train our LDA 
   MdlLinear = fitcdiscr(TrainData,TrainLabel);
   
   %test classify
   estLabel = predict(MdlLinear,TestData);
   
   %store classification
   stashClassification = [stashClassification; estLabel, TestLabel];
   stashData = [stashData; TestData];
end
fprintf('Done \n');

novstash = stashClassification(stashClassification(:,2) == 1,:);
corrnov = novstash(:,1) == novstash(:,2);
accnov = mean(corrnov)

expstash = stashClassification(stashClassification(:,2) == 3,:);
correxp = expstash(:,1) == expstash(:,2);
accexp = mean(correxp)


corr = stashClassification(:,1) == stashClassification(:,2);
acc = mean(corr)

% leave one out macro accuracies
% pegtx = 0.97;
% cutting = 0.94;
% suturing = 0.90;

% TaskAccuracy =
%     1.0000    0.9600    0.9231
%     0.8333    0.9000    0.8750

TaskAccuracy(:,tsk) = [accnov; accexp];
CountAccuracy(:,tsk) = [noviceCount; expertCount];

figure
gscatter3(Data(:,fc(1)),Data(:,fc(2)),Data(:,fc(3)),groupStr(Labels));
xlabel(key.c.All(fc(1)))
ylabel(key.c.All(fc(2)))
zlabel(key.c.All(fc(3)))
title('True Class')

figure
gscatter3(stashData(:,1),stashData(:,2),stashData(:,3),groupStr(stashClassification(:,1)));
xlabel(key.c.All(fc(1)))
ylabel(key.c.All(fc(2)))
zlabel(key.c.All(fc(3)))
title('Est Class')

