%% The Original EDGE study consisted of data collected on the Simulab EDGE platform
% at three different sites for three different FLS tasks.  Subjects were
% categorized for skill in different ways;
%  All data and results are combined into a matlab object via cells and
%  arrays and a human-readable indexing scheme. This object is:
% DataGlb =
%
%       dataLog: {1x583 cell}; The actual time-series data; each cell is an
%           individual run (iteration) of a single task by as single subject
%       content: {583x158 cell};  % contains all metadata info from the Excel file
%       Tasks: {'PegTx'  'Cutting'  'Suturing'}  % enumerates available
%         tasks
%       LogIdx: {[193x1 double]  [165x1 double]  [89x1 double]}; % the
%           indices  you use to acces dataLog{} cells
%       grp: [1x1 struct] % enumerates groups (e.g. skills) of logs and
%           partitions thereoff for k-fold cross-validation
%       grpFLS: [1x1 struct] % ignore this
%       lookupCol: % A function to help search .content by know key's (e.g.
%           video-filename.


%% For segment analysis, use the following:
%%%% How it works.... (internal vars)
% % % % Events are a change from one state to another: Open to Closed.  They can
% % % % be counted, e.g., sum(Events==Open); or serve to segment data: data(Events==Open)
% % % seg.GrL    = 1; % GraspVars (Fg,Qg) combined to extract a "grasp"; Open or Closed For all samples
% % % seg.GrEvL  = 2; % 'Events' assoc. with Gr; Samples: Open, Closed, or -1 (previous val)
% % % seg.GrR    = 3; % GraspVars (Fg,Qg) combined to extract a "grasp"; Open or Closed For all samples
% % % seg.GrEvR  = 4; % 'Events' assoc. with Gr; Samples: Open, Closed, or -1 (previous val)
% % % seg.FgL    = 5; % Only Fg-threshold based grasp detection: Open, Closed
% % % seg.FgEvL  = 6; % 'Events' assoc. with Fg; Open, Closed, -1(prev)
% % % seg.FgR    = 7; % Only Fg-threshold based grasp detection: Open, Closed
% % % seg.FgEvR  = 8; % 'Events' assoc. with Fg; Open, Closed, -1(prev)
% % % seg.ZSpdL  = 9; % zero-velocity segments: Stopped, Moving
% % % seg.ZSpdEvL=10; % 'Events' assoc with zVel: Stopped, Moving, -1(prev)
% % % seg.ZSpdR  =11; % zero-velocity segments: Stopped, Moving
% % % seg.ZSpdEvR=12; % 'Events' assoc with zVel: Stopped, Moving, -1(prev)
% % % 
% % % seg.Labels={...
% % %     'GrL', 'GrEvL', 'GrR', 'GrEvR',...
% % %     'FgL', 'FgEvL', 'FgR', 'FgEvR',...
% % %     'ZSpdL', 'ZSpdEvL', 'ZSpdR', 'ZSpdEvR'};
% % % 
% % % %seg.infoNumSeg   % total number of segments (cols in the info matrix)
% % % Seg.info.iN=1;    % st row of info matrix: number of samples (i) in this (col) segment
% % % Seg.info.iStart=2;% nd row of info matrix: index of start of this segment
% % % Seg.info.iEnd=3;  % rd row of info matrix: index of END of this segment
% % % Seg.info.pathLen=4;%th row of info matrix: path length of this segment
% % % Seg.info.disp  =5;% th row of info matrix: displacement (distance) beteen start and end point in cartesian space
% % % Seg.info.timeLen=6;%th row of info matrix: time duration of this segment
% % % Seg.info.EoM=7;   % th row of info matrix: EoM Economy of Motion (pathLen/Time) of this segment
% % % Seg.info.EoP=8;   % th row of info matrix: EoP " of Path (  disp /Time) of this segment
% % % Seg.info.EoD=9;   % th row of info matrix: EoD " of Displacement (  disp /Time) of this segment
% % % Seg.info.spdAtEnd=10;% th " "  info matrix: Mean speed from last N cm of segment
% % % %spdAtEndCmTh=2.0;% Threshold of last N cm to compute mean speed at end of segment for.
% % % 
% % % Seg.tag = fieldnames(Seg.info);

%% Segmentation Schemes (combine criteria of seg.* and seg.Labels to create
% % segmentation schemes segS.* like "between actuation events", "during
% % actuation (closed grasper)", or "between zero-speed segments" or
% % combinations therin.
% 
% SegS.betweenGrActs=1;   % All segments between Grasp Actuation Events (segs where grasper is open) (Based on Fg AND Qg)
% SegS.withinGrActs =2;   % All segments during GraspActuation Events (segs where grasper is closed)
% SegS.betweenFgActs=3;   % All segs between Fg-only-based grasp Actuation Events (segs when grasp force is low-open grasper)
% SegS.withinFgActs =4;   % All segs during " " " " (when grasp force is nonzero (above threshold): closed grasper on object)
% SegS.betweenZSpd  =5;   % All segs between zero-speed events
% SegS.Labels= {'betweenGrActs' 'withinGrActs' 'betweenFgActs' 'withinFgActs' 'betweenZSpd' };
% L=1;
% R=2;


%%
%load('D:\temp\EDGE_Segments.mat')

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

%FLS scores
flsColumn = DataGlb.lookupCol('FLS-score');
flsScores = cell2mat(DataGlb.content(DataGlb.LogIdx{tsk},flsColumn));
timeColumn =  DataGlb.lookupCol('General:TotalTime-homebase');
DataGlb.content(2,84)
DataGlb.content(3,84)

Hands = [1,2];
HandString = {'left', 'right'};


%To avoid having a ton of novice data
allsurgeons = [length(DataGlb.grp.all{tsk}.Idx{myGroups(1)}),length(DataGlb.grp.all{tsk}.Idx{myGroups(2)}),length(DataGlb.grp.all{tsk}.Idx{myGroups(3)})]
min_num_surgeons = min(allsurgeons)


%% Transform each segment to a new coordinate system: Origin is the end of the segment



%Clear this bad boy
SegData = [];
All = [];

tlength = [];
samplength = [];
totalSegs = [];

storeTime = [];
storeNumberSegments = [];

fprintf('Getting Segment Data...');
for gg = 1:length(myGroups) %looping through array of logs
    % Initialize
    tL = [];
    xleft = [];
    yleft = [];
    zleft = [];
    dxleft = [];
    dyleft = [];
    dzleft = [];
    ddxleft = [];
    ddyleft = [];
    ddzleft = [];
    dddxleft = [];
    dddyleft = [];
    dddzleft = [];
    gripL = [];
    dgripL = [];
    velmagL = [];
    velalphaL = [];
    velbetaL = [];
    accmagL = [];
    accalphaL = [];
    accbetaL = [];
    tR = [];
    xright = [];
    yright = [];
    zright = [];
    dxright = [];
    dyright = [];
    dzright = [];
    ddxright = [];
    ddyright = [];
    ddzright = [];
    dddxright = [];
    dddyright = [];
    dddzright = [];
    gripR = [];
    dgripR = [];
    velmagR = [];
    velalphaR = [];
    velbetaR = [];
    accmagR = [];
    accalphaR = [];
    accbetaR = [];
    
    %segment index
    segidx = 1;
    
    
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:allsurgeons(gg) %min_num_surgeons
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        %Get FLS score for that
        FLS_Trial = cell2mat(DataGlb.content(i,flsColumn));
        SegData{gg}.Trial{ii}.FLSScore = FLS_Trial;
        
        
        % ADD SEGMENTS:
        % Extract actual segment data: (left first)
        mySegmentsL = SegScheme{tsk}.dataLogSegs(i, SegType ,L); % Left hand
        
        for s = 1:length(mySegmentsL.sgidx) %loop through all segments for that particular task   
            
            % selects only the samples of this segment and get start/end
            iSamples = mySegmentsL.sgidx{s};
            iEnd = iSamples(end);
            iStart = iSamples(1);
            
            %For each segment
            %position
            t_seg = DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iStart, G.Time);
            xleft_seg = DataGlb.dataLog{i}( iSamples, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL);
            yleft_seg = DataGlb.dataLog{i}( iSamples, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL);
            zleft_seg = DataGlb.dataLog{i}( iSamples, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL);

            %velocity
            dxleft_seg = DataGlb.dataLog{i}( iSamples, G.dxL) ;
            dyleft_seg = DataGlb.dataLog{i}( iSamples, G.dyL) ;
            dzleft_seg = DataGlb.dataLog{i}( iSamples, G.dzL) ;

            %accel
            ddxleft_seg = DataGlb.dataLog{i}( iSamples, G.ddxL) ;
            ddyleft_seg = DataGlb.dataLog{i}( iSamples, G.ddyL) ;
            ddzleft_seg = DataGlb.dataLog{i}( iSamples, G.ddzL) ;
            %jerk
            dddxleft_seg = DataGlb.dataLog{i}( iSamples, G.dddxL) ;
            dddyleft_seg = DataGlb.dataLog{i}( iSamples, G.dddyL) ;
            dddzleft_seg = DataGlb.dataLog{i}( iSamples, G.dddzL) ;
            
            %grp pos and vel
            tempGL = DataGlb.dataLog{i}(iSamples,G.QgL);
            tempdGL = DataGlb.dataLog{i}(iSamples,G.dQgL);
            
            %angles
            %magnitude and angle of velocity
            tmpvelmag = NormRowWise([dxleft_seg,dyleft_seg,dzleft_seg]);
            [tmpvelalpha,tmpvelbeta] = PitchYaw3D([dxleft_seg,dyleft_seg,dzleft_seg]);

            %magnitude and angle of acceleration
            tmpaccmag = NormRowWise([ddxleft_seg,ddyleft_seg,ddzleft_seg]);
            [tmpaccalpha,tmpaccbeta] = PitchYaw3D([ddxleft_seg,ddyleft_seg,ddzleft_seg]);

            %add to grouping
            tL = [ tL ; t_seg] ;
            gripL = [gripL; tempGL];
            dgripL = [dgripL; tempdGL];
            xleft = [ xleft ; xleft_seg] ;
            yleft = [ yleft ; yleft_seg] ;
            zleft = [ zleft ; zleft_seg] ;

            dxleft = [ dxleft ; dxleft_seg ] ;
            dyleft = [ dyleft ; dyleft_seg ] ;
            dzleft = [ dzleft ; dzleft_seg ] ;
            
            ddxleft = [ddxleft ; ddxleft_seg];
            ddyleft = [ddyleft ; ddyleft_seg];
            ddzleft = [ddzleft ; ddzleft_seg];
            
            dddxleft = [dddxleft ; dddxleft_seg];
            dddyleft = [dddyleft ; dddyleft_seg];
            dddzleft = [dddzleft ; dddzleft_seg];
            
            velmagL = [velmagL; tmpvelmag];
            velalphaL = [velalphaL; tmpvelalpha];
            velbetaL = [velbetaL; tmpvelbeta];
            
            accmagL = [accmagL; tmpaccmag];
            accalphaL = [accalphaL; tmpaccalpha];
            accbetaL = [accbetaL; tmpaccbeta];

            %store individual segments
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Time = t_seg;
            
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Pos = [xleft_seg, yleft_seg, zleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Vel = [dxleft_seg, dyleft_seg, dzleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Acc = [ddxleft_seg, ddyleft_seg, ddzleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Jerk = [dddxleft_seg, dddyleft_seg, dddzleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Grip = [tempGL ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.dGrip = [tempdGL ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.VelAngle = [tmpvelmag, tmpvelalpha, tmpvelbeta];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.AccAngle = [tmpaccmag, tmpaccalpha, tmpaccbeta ];
        
            segidx = segidx +1;
        end
        
        mySegmentsR = SegScheme{tsk}.dataLogSegs(i, SegType ,R); % Now Right hand
        
        for s = 1:length(mySegmentsR.sgidx) %loop through all segments for that particular task   
            
            % selects only the samples of this segment and get start/end
            iSamples = mySegmentsR.sgidx{s};
            iEnd = iSamples(end);
            iStart = iSamples(1);
            
            %For each segment
            t_seg = DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iStart, G.Time);
            %position
            xright_seg = DataGlb.dataLog{i}( iSamples, G.xR) - DataGlb.dataLog{i}( iEnd, G.xR);
            yright_seg = DataGlb.dataLog{i}( iSamples, G.yR) - DataGlb.dataLog{i}( iEnd, G.yR);
            zright_seg = DataGlb.dataLog{i}( iSamples, G.zR) - DataGlb.dataLog{i}( iEnd, G.zR);
            
            %velocity
            dxright_seg = DataGlb.dataLog{i}( iSamples, G.dxR) ;
            dyright_seg = DataGlb.dataLog{i}( iSamples, G.dyR) ;
            dzright_seg = DataGlb.dataLog{i}( iSamples, G.dzR) ;

            %acceleration
            ddxright_seg = DataGlb.dataLog{i}( iSamples, G.ddxR) ;
            ddyright_seg = DataGlb.dataLog{i}( iSamples, G.ddyR) ;
            ddzright_seg = DataGlb.dataLog{i}( iSamples, G.ddzR) ;
            
            %jerk
            dddxright_seg = DataGlb.dataLog{i}( iSamples, G.dddxR) ;
            dddyright_seg = DataGlb.dataLog{i}( iSamples, G.dddyR) ;
            dddzright_seg = DataGlb.dataLog{i}( iSamples, G.dddzR) ;
            
            %angles
            %magnitude and angle of velocity
            tmpvelmag = NormRowWise([dxright_seg,dyright_seg,dzright_seg]);
            [tmpvelalpha,tmpvelbeta] = PitchYaw3D([dxright_seg,dyright_seg,dzright_seg]);

            %magnitude and angle of acceleration
            tmpaccmag = NormRowWise([ddxright_seg,ddyright_seg,ddzright_seg]);
            [tmpaccalpha,tmpaccbeta] = PitchYaw3D([ddxright_seg,ddyright_seg,ddzright_seg]);
            
            %gripper angle and vel
            tempGR = DataGlb.dataLog{i}(iSamples,G.QgR);
            tempdGR = DataGlb.dataLog{i}(iSamples,G.dQgR);
            
            %add to grouping
            tR = [ tR ; t_seg] ;
            gripR = [gripR; tempGR];
            dgripR = [dgripR; tempdGR];
            xright = [ xright ; xright_seg] ;
            yright = [ yright ; yright_seg] ;
            zright = [ zright ; zright_seg] ;

            dxright = [ dxright ; dxright_seg ] ;
            dyright = [ dyright ; dyright_seg ] ;
            dzright = [ dzright ; dzright_seg ] ;
            
            ddxright = [ddxright ; ddxright_seg];
            ddyright = [ddyright ; ddyright_seg];
            ddzright = [ddzright ; ddzright_seg];
            
            dddxright = [dddxright ; dddxright_seg];
            dddyright = [dddyright ; dddyright_seg];
            dddzright = [dddzright ; dddzright_seg];
            
            velmagR = [velmagR; tmpvelmag];
            velalphaR = [velalphaR; tmpvelalpha];
            velbetaR = [velbetaR; tmpvelbeta];
            
            accmagR = [accmagR; tmpaccmag];
            accalphaR = [accalphaR; tmpaccalpha];
            accbetaR = [accbetaR; tmpaccbeta];


            %store individual segments
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Time = t_seg;
            
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Pos = [xright_seg, yright_seg, zright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Vel = [dxright_seg, dyright_seg, dzright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Acc = [ddxright_seg, ddyright_seg, ddzright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Jerk = [dddxright_seg, dddyright_seg, dddzright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Grip = [tempGR ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.dGrip = [tempdGR ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.VelAngle = [tmpvelmag, tmpvelalpha, tmpvelbeta];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.AccAngle = [tmpaccmag, tmpaccalpha, tmpaccbeta ];
        
            segidx = segidx +1;
        end
        
        
        storeTime = [storeTime; cell2mat(DataGlb.content(i,timeColumn))];
        %stash the number of Left hand and Right Hand segments
        storeNumberSegments = [storeNumberSegments; length(mySegmentsL.sgidx), length(mySegmentsR.sgidx)];
    end
    
    %Store all
    PosL = [xleft, yleft, zleft ];
    VelL = [dxleft, dyleft, dzleft ];
    AccL = [ddxleft, ddyleft, ddzleft ];
    JerkL = [dddxleft, dddyleft, dddzleft ];
    
    PosR = [xright, yright, zright ];
    VelR = [dxright, dyright, dzright ];
    AccR = [ddxright, ddyright, ddzright ];
    JerkR = [dddxright, dddyright, dddzright ];
    
    All{gg}.TimeL = tL;
    All{gg}.VelL = VelL;
    All{gg}.AccL = AccL;
    All{gg}.JerkL = JerkL;
    All{gg}.PosL = PosL;
    All{gg}.GripL = gripL ;
    All{gg}.dGripL = dgripL ;
    All{gg}.VelAngleL = [velmagL, velalphaL, velbetaL];
    All{gg}.AccAngleL = [accmagL, accalphaL, accbetaL ];

    All{gg}.TimeR = tR;
    All{gg}.VelR = VelR;
    All{gg}.AccR = AccR;
    All{gg}.JerkR = JerkR;
    All{gg}.PosR = PosR;
    All{gg}.GripR = gripR ;
    All{gg}.dGripR = dgripR ;
    All{gg}.VelAngleR = [velmagR, velalphaR, velbetaR];
    All{gg}.AccAngleR = [accmagR, accalphaR, accbetaR ];

    %store how many total segments we have
    totalSegs(gg) = segidx;
end


mean(storeTime)
sum(mean(storeNumberSegments))
mean(std(storeNumberSegments))

totalSegs

fprintf('Done \n');

%% Create Data and labels

%key for later
key = [];
key.col.dx = 1;
key.col.dy = 2;
key.col.dz = 3;
key.col.ddx = 4;
key.col.ddy = 5;
key.col.ddz = 6;
key.col.dddx = 7;
key.col.dddy = 8;
key.col.dddz = 9;
key.col.velmag = 10;
key.col.velalpha = 11;
key.col.velbeta = 12;
key.col.accmag = 13;
key.col.accalpha = 14;
key.col.accbeta = 15;
key.col.grip = 16;
key.col.dgrip = 17;

key.strings = {'dx','dy','dz','ddx','ddy','ddz','dddx','dddy','dddz','velmag','velalpha','velbeta','accmag','accalpha','accbeta','grip','dgrip'};

%which states we want to keep
StrLabels = {'Novice','Expert'};
grplist = [1,3];

%clear this monster
DataTrain = [];
DataValidate = [];
DataTest = [];
LabelsTrain = [];
LabelsValidate = [];
LabelsTest = [];

%sample into training/validation/testing
trainRatio = 0.7;
valRatio = 0.2;
testRatio = 0.1;

%segment data
AllSegmentCount = [];
AllSegmentDuration = [];

fprintf('Storing into a matrix for training...');
for gg = grplist %skill levels
    %for sub sampling Training, Validation and Testing
    [idxtrain, idxval, idxtest] = dividerand(totalSegs(gg),trainRatio,valRatio,testRatio);
    segmentid = 1;
    
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        segz = length(SegData{gg}.Trial{ii}.Hand{1}.Segment) + length(SegData{gg}.Trial{ii}.Hand{2}.Segment);
        AllSegmentCount =[AllSegmentCount; segz, gg];
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            
            
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};

                AllSegmentDuration = [AllSegmentDuration; length(struct.Time), gg];

                %stash it all in a temporary mat
                dtemp = [];
                dtemp(:,key.col.dx) = struct.Vel(:,1);
                dtemp(:,key.col.dy) = struct.Vel(:,2);
                dtemp(:,key.col.dz) = struct.Vel(:,3);
                dtemp(:,key.col.ddx) = struct.Acc(:,1);
                dtemp(:,key.col.ddy) = struct.Acc(:,2);
                dtemp(:,key.col.ddz) = struct.Acc(:,3);
                dtemp(:,key.col.dddx) = struct.Jerk(:,1);
                dtemp(:,key.col.dddy) = struct.Jerk(:,2);
                dtemp(:,key.col.dddz) = struct.Jerk(:,3);
                dtemp(:,key.col.grip) = struct.Grip;
                dtemp(:,key.col.dgrip) = struct.dGrip;
                
                dtemp(:,key.col.velmag) = struct.VelAngle(:,1);
                dtemp(:,key.col.velalpha) = struct.VelAngle(:,2);
                dtemp(:,key.col.velbeta) = struct.VelAngle(:,3);
                dtemp(:,key.col.accmag) = struct.AccAngle(:,1);
                dtemp(:,key.col.accalpha) = struct.AccAngle(:,2);
                dtemp(:,key.col.accbeta) = struct.AccAngle(:,3);
                
                ng = size(dtemp,1);

                if(any(segmentid==idxtrain))
                    DataTrain = [DataTrain; dtemp];
                    LabelsTrain = [LabelsTrain; ones(ng,1)*gg];
                elseif(any(segmentid==idxval))
                    DataValidate = [DataValidate; dtemp];
                    LabelsValidate = [LabelsValidate; ones(ng,1)*gg];
                elseif(any(segmentid==idxtest))
                    DataTest = [DataTest; dtemp];
                    LabelsTest = [LabelsTest; ones(ng,1)*gg];
                end
                segmentid = segmentid + 1;
            end
        end
    end
end
fprintf('Done \n');

%test columns
tc1 = key.col.dx;
tc2 = key.col.dy;
tc3 = key.col.dz;
% tc1 = key.col.velmag;
% tc2 = key.col.velalpha;
% tc3 = key.col.velbeta;
% tc1 = key.col.accmag;
% tc2 = key.col.accalpha;
% tc3 = key.col.accbeta;

if(false)
figure
gscatter3(DataTrain(:,tc1),DataTrain(:,tc2),DataTrain(:,tc3),LabelsTrain)
title('phase portrait train')

figure
gscatter3(DataTest(:,tc1),DataTest(:,tc2),DataTest(:,tc3),LabelsTest)
title('phase portrait test')
end

%check our ratios
allsizes = [size(DataTrain,1),size(DataValidate,1),size(DataTest,1)];
sumsizes = sum(allsizes);
actualratios = allsizes ./ sumsizes

novSegmentDuration = AllSegmentDuration(AllSegmentDuration(:,2) == 1,1);
expSegmentDuration = AllSegmentDuration(AllSegmentDuration(:,2) == 3,1);

mean(novSegmentDuration)
mean(expSegmentDuration)

novSegmentCount = AllSegmentCount(AllSegmentCount(:,2) == 1,1);
expSegmentCount = AllSegmentCount(AllSegmentCount(:,2) == 3,1);

mean(novSegmentCount)
std(novSegmentCount)
mean(expSegmentCount)
std(expSegmentCount)

%% SCALE ALL DIMENSIONS TO HAVE THE SAME RANGE (See Lin2006 Paper)

DataTrainRaw = DataTrain;
DataTestRaw = DataTest;
DataValidateRaw = DataValidate;

%Scaling function based on variance and mean
DataTrain = NormalizeFeatures(DataTrainRaw);
DataTest = NormalizeFeatures(DataTestRaw);
DataValidate = NormalizeFeatures(DataValidateRaw);

%test plot
plotcols2 = [key.col.ddx, key.col.velmag, key.col.accmag]

figure
gscatter3(DataTrainRaw(:,tc1), DataTrainRaw(:,tc2), DataTrainRaw(:,tc3),LabelsTrain)
xlabel(key.strings(plotcols2(1)));
ylabel(key.strings(plotcols2(2)));
zlabel(key.strings(plotcols2(3)));
title('UnScaled Feature Data')

figure
gscatter3(DataTrain(:,tc1), DataTrain(:,tc2), DataTrain(:,tc3),LabelsTrain)
xlabel(key.strings(plotcols2(1)));
ylabel(key.strings(plotcols2(2)));
zlabel(key.strings(plotcols2(3)));
title('Scaled Feature Data')

%%
%save the data and labels to a mat for easy loading

save DataMatrixEdge.mat DataTrain LabelsTrain DataValidate DataTrainRaw DataTestRaw DataValidateRaw LabelsValidate DataTest LabelsTest key
%load DataMatrixEdge.mat

%% Run relief f on all states cuz why the hell not

[RANK,WEIGHT] = relieff(DataTrain,LabelsTrain,100);
fsize = 14
RANK
WEIGHT(RANK)
beststates = key.strings(RANK)
%These are the new states with scaled states
%RANK = [16,9,8,17,10,7,13,3,2,1,5,4,6,11,12,14,15]
%strings = {'grip','dddz','dddy','dgrip','velmag','dddx','accmag','dz','dy','dx','ddy','ddx','ddz','velalpha','velbeta','accalpha','accbeta'}
%weights = [0.00233821624040019,0.00277628860437168,0.00301852920541321,0.00198931956255718,0.00227216502623299,0.00135308611103428,0.00484540834484069,0.00567846559114360,0.00686884662561775,0.00496121156556170,0.00110351075253696,-0.000114448607495619,0.00378187239168139,-0.000963230401998722,-0.00211511503819641,0.00759164985191833,0.00544718082661790]

%%

% plotcols = [RANK(1),RANK(2),RANK(3)];
plotcols = [16,9,8];
beststates = key.strings(plotcols)
strgrp = {'Novice';'Intermeditate';'Expert'}
figure
h = gscatter3(DataTrainRaw(:,plotcols(1)),DataTrainRaw(:,plotcols(2)),DataTrainRaw(:,plotcols(3)),strgrp(LabelsTrain)',{'c','r'},{'o','+'}); %
xlabel('$\theta$','FontSize',fsize,'interpreter','latex');
% ylabel('$\frac{d^{3}Z}{dt^{3}}$','FontSize',fsize,'interpreter','latex');
ylabel('$\ddot{Z}$','FontSize',fsize,'interpreter','latex')
% zlabel('$\frac{d^{3}Y}{dt^{3}}$','FontSize',fsize,'interpreter','latex');
zlabel('$\ddot{Y}$','FontSize',fsize,'interpreter','latex')
%title('Top RELIEFF states','FontSize',fsize)
leg = findobj(gcf,'Tag','legend');
set(leg,'FontSize',fsize);


%% Find best states for Seperability (RUN ME AT SCHOOL)

testcolumns = [key.col.grip,key.col.dgrip,key.col.ddx,key.col.ddy,key.col.ddz,key.col.velmag,key.col.dddx,key.col.dddy,key.col.dddz,key.col.accmag];
%testcolumns = [key.col.velmag, key.col.velalpha, key.col.velbeta, key.col.accmag, key.col.accalpha, key.col.accbeta, key.col.ddx, key.col.ddz];
%testcolumns = [key.col.grip,key.col.dgrip,key.col.dddx,key.col.velmag,key.col.accmag,key.col.dddz,key.col.dddy,key.col.ddx];
nkcombos = nchoosekmulti(testcolumns,2:length(testcolumns))


[Variations,BestSep,BestSepClass,Data_All,Diff_All] = RelativeRBFSeperabilityCheck(DataTrain,LabelsTrain,testcolumns,key.strings,1,'ploton');
BestSep.bestcollabels
BestSep.max
idn = 3;
BestSepClass{idn}.max
BestSepClass{idn}.metric
columnrbf = BestSepClass{idn}.bestcolumns
collabs = BestSepClass{idn}.bestcollabels

%these are the actual best
columnrbf = [16,7,9];
beststates = key.strings(columnrbf)

%THESE ARE THE NEW SCALED STATES
%for 4D, best = ['grip'    'dddx'    'dddy'    'dddz'] = [14,7,8,9]
%with weight = 0.0058

%for 3D, best = ['grip'    'dddx'    'dddz'] = [16,7,9]
%with weight = 0.0063

%for 2D, best = ['grip'    'dddz'] = [16,9]
%with weight = 0.0068

%for 1D, best = ['grip'] = [16]
%with weight = 0.0097

if(length(columnrbf) < 3)
    
    figure
    gscatter(DataTrain(:,columnrbf(1)),DataTrain(:,columnrbf(2)),LabelsTrain);
    title('class data')
     
    figure
    gscatter(DataTrain(:,columnrbf(1)),DataTrain(:,columnrbf(2)),LabelsTrain);
    hold on
    cc = 1;
    Surface3D(BestSep.DataClass{cc}(:,1),BestSep.DataClass{cc}(:,2),BestSep.Difference{cc});
    hold on
    cc = 2;
    Surface3D(BestSep.DataClass{cc}(:,1),BestSep.DataClass{cc}(:,2),BestSep.Difference{cc});
    hold off
    title('best seperation')
    colormap cool
    colorbar
    xlabel(key.strings(columnrbf(1)));
    ylabel(key.strings(columnrbf(2)));
    zlabel('seperability');
end

%%
columnrbf = [16,7,9];
beststates = key.strings(columnrbf)
strgrp = {'Novice';'Intermeditate';'Expert'}

figure
h = gscatter3(DataTrainRaw(:,columnrbf(1)),DataTrainRaw(:,columnrbf(2)),DataTrainRaw(:,columnrbf(3)),strgrp(LabelsTrain)',{'c','r'},{'o','+'});

xlabel('$\theta$','FontSize',fsize,'interpreter','latex');
%\dddot doesnt exist
% ylabel('$\frac{d^{3}X}{dt^{3}}$','FontSize',fsize,'interpreter','latex');
ylabel('$\ddot{X}$','FontSize',fsize,'interpreter','latex')
% zlabel('$\frac{d^{3}Z}{dt^{3}}$','FontSize',fsize,'interpreter','latex');
zlabel('$\ddot{Z}$','FontSize',fsize,'interpreter','latex')
%title('Top RELIEFRBF states','FontSize',fsize)
leg = findobj(gcf,'Tag','legend');
set(leg,'FontSize',fsize);
    
%% sep as color plots

figure
cc = 1;
scatter3(BestSepClass{idn}.DataClass{cc}(:,1),BestSepClass{idn}.DataClass{cc}(:,2),BestSepClass{idn}.DataClass{cc}(:,3),10,BestSepClass{idn}.Difference{cc});
hold on
cc = 2;
scatter3(BestSepClass{idn}.DataClass{cc}(:,1),BestSepClass{idn}.DataClass{cc}(:,2),BestSepClass{idn}.DataClass{cc}(:,3),10,BestSepClass{idn}.Difference{cc});
hold off
xlabel(key.strings(columnrbf(1)));
ylabel(key.strings(columnrbf(2)));
zlabel(key.strings(columnrbf(3)));
title('seperability as color')
colormap cool
%colorbar
hc = colorbar;
ylabel(hc, 'W_{rbf}','FontSize',fsize)


%% lets try a random forest with these best labels

%Results Accuracy 70.5% with leave surgeon out
%Out of bag error: 0.2742

%test columns
testcols = columnrbf;
%this is for
%best = ['grip'    'dddx'    'dddz'] = [16,7,9]
%with weight = 0.0063

rng(1); % For reproducibility
Mdl = TreeBagger(100,DataTrain(:,testcols),LabelsTrain,'OOBPrediction','On',...
    'Method','classification')

% view(Mdl.Trees{1},'Mode','graph') %view it but its huge

%plot loss function
figure;
oobErrorBaggedEnsemble = oobError(Mdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';

labelestimate = predict(Mdl,DataTest(:,testcols));

corr = str2num(cell2mat(labelestimate)) == LabelsTest;

acc = mean(corr)

oobmin = min(oobErrorBaggedEnsemble)

return;

%% compare difference between velalpha and velbeta ITS BAD

DataPlot = DataTrain(:,[key.col.velalpha,key.col.velbeta]);

DiffABvel = abs(DataTrain(:,key.col.velalpha) - DataTrain(:,key.col.velbeta));
DiffABacc = abs(DataTrain(:,key.col.accalpha) - DataTrain(:,key.col.accbeta));

figure
gscatter(DiffABvel,DiffABacc,LabelsTrain);
title('diff in angles')

figure
gscatter(DataPlot(:,1),DataPlot(:,2),LabelsTrain);
title('class data')

%% Test RANDOM FORESTS YAY

%test columns
% testcols = [key.col.velmag, key.col.velalpha, key.col.velbeta, key.col.accmag, key.col.accalpha, key.col.accbeta];
testcols = [key.col.velmag, key.col.accmag];


rng(1); % For reproducibility
Mdl = TreeBagger(50,DataTrain(:,testcols),LabelsTrain,'OOBPrediction','On',...
    'Method','classification')

% view(Mdl.Trees{1},'Mode','graph') %view it but its huge

%plot loss function
figure;
oobErrorBaggedEnsemble = oobError(Mdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';

labelestimate = predict(Mdl,DataTest(:,testcols));

corr = str2num(cell2mat(labelestimate)) == LabelsTest;

acc = mean(corr)

%% plot bounds of random forest

bounds = DataBounds(DataTest(:,testcols));
grid = ndimgrid(bounds,100);
labelg = predict(Mdl,grid);

figure
gscatter(grid(:,1),grid(:,2),labelg);
xlabel(key.strings(testcols(1)))
ylabel(key.strings(testcols(2)))
title('grid bounds')

figure
gscatter(DataTest(:,testcols(1)),DataTest(:,testcols(2)),labelestimate);
xlabel(key.strings(testcols(1)))
ylabel(key.strings(testcols(2)))
title('actual data')


%% simple rbf

%test columns
%testcols = [key.col.velmag, key.col.velalpha, key.col.velbeta];
testcols = [key.col.dx, key.col.dy, key.col.dz, key.col.ddx, key.col.ddy, key.col.ddz, key.col.velmag, key.col.velalpha, key.col.velbeta, key.col.accmag, key.col.accalpha, key.col.accbeta];

DataOn = DataTest(:,testcols);

[Difference,ClassData,ProbData] = SimpleRelativeRBFTrain(DataOn,LabelsTest);

figure
gscatter3(DataOn(:,1),DataOn(:,2),DataOn(:,3),LabelsTest);
title('classes')

figure
cc = 1;
scatter3(ClassData{cc}(:,1),ClassData{cc}(:,2),ClassData{cc}(:,3),10,Difference{cc});
hold on
cc = 2;
scatter3(ClassData{cc}(:,1),ClassData{cc}(:,2),ClassData{cc}(:,3),10,Difference{cc});
hold off
title('seperability')
colormap cool
colorbar




