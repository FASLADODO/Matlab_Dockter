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

tsk = 1; %PegTx:1 Cutting:2 Suturing:3
TaskStr = DataGlb.Tasks;

myGroups = [g.flsNov g.flsInt g.gtExp];
groupStr = {'Novice','Intermediate','Expert'};

% Choose Segment Type
%SegType = SegS.betweenFgActs; % (based only on Forces, excludes position
%SegType = SegS.betweenZSpd; % (once the tool slows to zero speed)
SegType = SegS.betweenGrActs; % (includes forces and grasper position)

Hands = [1,2];
HandString = {'left', 'right'};

%To avoid having a ton of novice data
allsurgeons = [length(DataGlb.grp.all{tsk}.Idx{myGroups(1)}),length(DataGlb.grp.all{tsk}.Idx{myGroups(2)}),length(DataGlb.grp.all{tsk}.Idx{myGroups(3)})]
min_num_surgeons = min(allsurgeons)


%% Access FLS score

%You can find the column headings for each of the columns in the second row of this struct
%Get column names
contentcolumns = DataGlb.content(2,:);% look at this to figure out string

%FLS scores
flsColumn = DataGlb.lookupCol('FLS-score');

%lap cases
casesColumn = DataGlb.lookupCol('AllTimeTotalLaprProceduresEstimate');
trainingColumn = DataGlb.lookupCol('TrainingLevel');

%task time
timeColumn = DataGlb.lookupCol('General:TotalTime-homebase');

% DataGlb.content{ix, type_of_data_you_want}

AllData = [];
for task = 1:3
    flsScores = cell2mat(DataGlb.content(DataGlb.LogIdx{tsk},flsColumn));
    cases = cell2mat(DataGlb.content(DataGlb.LogIdx{tsk},casesColumn));
    training = cell2mat(DataGlb.content(DataGlb.LogIdx{tsk},trainingColumn));
    time = cell2mat(DataGlb.content(DataGlb.LogIdx{tsk},timeColumn));
    AllData = [AllData; flsScores, cases, training, time];
end





