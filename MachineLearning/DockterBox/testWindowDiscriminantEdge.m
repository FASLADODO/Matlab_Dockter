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
    fprintf('Loading Global Data Structure EdgeDataGlb.mat (huge) ...');
    try
        load('C:\temp\EdgeDataGlb.mat')
    catch
        disp('Could not load local. Loading from M drive...');
        load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EdgeDataGlb.mat')
    end
    try
        load('C:\temp\EDGE_Segments.mat')
    catch
        disp('Could not load local. Loading from M drive...');
        load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EDGE_Segments.mat')
    end
    fprintf(' DONE.\n');
end



% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boundtype = 'hull'; %'hull'
rbfstates = 'velacc'; %'velacc'; 'accjerk'; 'velaccjerk'; 'posvelacc'

tsk = 1; %PegTx:1 Cutting:2 Suturing:3
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

ggLbl = {'r', 'g', 'b'};

Hands = [1,2];
HandString = {'left', 'right'};




% DEMO CODE %%%%%%%%%%%
% Chose Log.
i = DataGlb.LogIdx{tsk}(1);% my specific LogIdx (get this from grp scructure for experts, novices, etc and for loop it
% Extract Info about samples: (start and Stop samples)
iStartL = SegScheme{tsk}.info(i, SegType , Seg.info.iStart, L).d; % Left Hand
iStartR = SegScheme{tsk}.info(i, SegType , Seg.info.iStart, R).d; % Right Hand
iEndL = SegScheme{tsk}.info(i, SegType , Seg.info.iEnd, L).d; % Left Hand
iEndR = SegScheme{tsk}.info(i, SegType , Seg.info.iEnd, R).d; % Right Hand
... Select whatever info you want via Seg.info.<> struct
    
% Extract actual segment data:
mySegmentsL = SegScheme{tsk}.dataLogSegs(i, SegType ,L); % Left hand
mySegmentsR = SegScheme{tsk}.dataLogSegs(i, SegType ,R); % Right hand

figure
%access like this...get the segment index list
for s = 1:length(mySegmentsL.sgidx)
    
    % selects only the samples of this segment
    iSamples = mySegmentsL.sgidx{s}; 
    % do stuff with only those samples
    plot( DataGlb.dataLog{i}( iSamples, G.Time), ...
          DataGlb.dataLog{i}( iSamples, G.toolPathR),'.b'); hold on;

end





%% Transform each segment to a new coordinate system: Origin is the end of the segment


%Clear this bad boy
SegData = [];
All = [];
All{1}.Input = [];
All{1}.Speed = [];
All{1}.Acc = [];
All{1}.Jerk = [];
All{1}.Pos  = [];
All{1}.Grp = [];

All{2}.Input = [];
All{2}.Speed = [];
All{2}.Acc = [];
All{2}.Jerk = [];
All{2}.Pos  = [];
All{2}.Grp = [];

All{3}.Input = [];
All{3}.Speed = [];
All{3}.Acc = [];
All{3}.Jerk = [];
All{3}.Pos  = [];
All{3}.Grp = [];

tlength = [];
samplength = [];

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
    grpL = [];
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
    grpR = [];
    
    tskcnt = 0;
    id = 1;
    
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)})
        
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
            t_seg = DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iStart, G.Time);
            xleft_seg = DataGlb.dataLog{i}( iSamples, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL);
            yleft_seg = DataGlb.dataLog{i}( iSamples, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL);
            zleft_seg = DataGlb.dataLog{i}( iSamples, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL);

            dxleft_seg = DataGlb.dataLog{i}( iSamples, G.dxL) ;
            dyleft_seg = DataGlb.dataLog{i}( iSamples, G.dyL) ;
            dzleft_seg = DataGlb.dataLog{i}( iSamples, G.dzL) ;

            grpL_seg = repmat( DataGlb.grp.tag(myGroups(gg)), [length(iSamples),1]) ;
            grp_num = repmat( gg, [length(iSamples),1]) ;
            
            T = mean(diff(t_seg)); %for computing accelerations

            %ddxleft_seg = Calculate_velocity( dxleft_seg, T, 'holobrodko') ;
            %ddyleft_seg = Calculate_velocity( dyleft_seg, T, 'holobrodko') ;
            %ddzleft_seg = Calculate_velocity( dzleft_seg, T, 'holobrodko') ;
            %dddxleft_seg = Calculate_velocity( ddxleft_seg, T, 'holobrodko') ;
            %dddyleft_seg = Calculate_velocity( ddyleft_seg, T, 'holobrodko') ;
            %dddzleft_seg = Calculate_velocity( ddzleft_seg, T, 'holobrodko') ;
            
            ddxleft_seg = DataGlb.dataLog{i}( iSamples, G.ddxL) ;
            ddyleft_seg = DataGlb.dataLog{i}( iSamples, G.ddyL) ;
            ddzleft_seg = DataGlb.dataLog{i}( iSamples, G.ddzL) ;
            dddxleft_seg = DataGlb.dataLog{i}( iSamples, G.dddxL) ;
            dddyleft_seg = DataGlb.dataLog{i}( iSamples, G.dddyL) ;
            dddzleft_seg = DataGlb.dataLog{i}( iSamples, G.dddzL) ;
            
            tlength{gg}(id) = t_seg(end) - t_seg(1);
            samplength{gg}(id) = length(iSamples);
            id = id +1;
            
            
            
            %add to grouping
            tL = [ tL ; t_seg] ;
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

            grpL = [ grpL; grpL_seg ];

            %store individual segments
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Time = t_seg;
            
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Pos = [xleft_seg, yleft_seg, zleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Vel = [dxleft_seg, dyleft_seg, dzleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Acc = [ddxleft_seg, ddyleft_seg, ddzleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Jerk = [dddxleft_seg, dddyleft_seg, dddzleft_seg ];
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.RBF = [];
            
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.Grp = grpL_seg;
            SegData{gg}.Trial{ii}.Hand{L}.Segment{s}.GrpNum = grp_num;
            
            tskcnt = tskcnt + 1;
        end
        
        mySegmentsR = SegScheme{tsk}.dataLogSegs(i, SegType ,R); % Now Right hand
        
        for s = 1:length(mySegmentsR.sgidx) %loop through all segments for that particular task   
            
            % selects only the samples of this segment and get start/end
            iSamples = mySegmentsR.sgidx{s};
            iEnd = iSamples(end);
            iStart = iSamples(1);
            
            %For each segment
            t_seg = DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iStart, G.Time);
            xright_seg = DataGlb.dataLog{i}( iSamples, G.xR) - DataGlb.dataLog{i}( iEnd, G.xR);
            yright_seg = DataGlb.dataLog{i}( iSamples, G.yR) - DataGlb.dataLog{i}( iEnd, G.yR);
            zright_seg = DataGlb.dataLog{i}( iSamples, G.zR) - DataGlb.dataLog{i}( iEnd, G.zR);

            dxright_seg = DataGlb.dataLog{i}( iSamples, G.dxR) ;
            dyright_seg = DataGlb.dataLog{i}( iSamples, G.dyR) ;
            dzright_seg = DataGlb.dataLog{i}( iSamples, G.dzR) ;

            grpR_seg = repmat( DataGlb.grp.tag(myGroups(gg)), [length(iSamples),1]) ;
            grp_num = repmat( gg, [length(iSamples),1]) ;
            
            T = mean(diff(t_seg)); %for computing accelerations
            %ddxright_seg = Calculate_velocity( dxright_seg, T, 'holobrodko') ;
            %ddyright_seg = Calculate_velocity( dyright_seg, T, 'holobrodko') ;
            %ddzright_seg = Calculate_velocity( dzright_seg, T, 'holobrodko') ;
            %dddxright_seg = Calculate_velocity( ddxright_seg, T, 'holobrodko') ;
            %dddyright_seg = Calculate_velocity( ddyright_seg, T, 'holobrodko') ;
            %dddzright_seg = Calculate_velocity( ddzright_seg, T, 'holobrodko') ;
            
            ddxright_seg = DataGlb.dataLog{i}( iSamples, G.ddxR) ;
            ddyright_seg = DataGlb.dataLog{i}( iSamples, G.ddyR) ;
            ddzright_seg = DataGlb.dataLog{i}( iSamples, G.ddzR) ;
            dddxright_seg = DataGlb.dataLog{i}( iSamples, G.dddxR) ;
            dddyright_seg = DataGlb.dataLog{i}( iSamples, G.dddyR) ;
            dddzright_seg = DataGlb.dataLog{i}( iSamples, G.dddzR) ;
            
            tlength{gg}(id) = t_seg(end) - t_seg(1);
            samplength{gg}(id) = length(iSamples);
            id = id +1;
             
            

            %add to grouping
            tR = [ tR ; t_seg] ;
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

            grpR = [ grpR; grpR_seg ];

            %store individual segments
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Time = t_seg;
            
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Pos = [xright_seg, yright_seg, zright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Vel = [dxright_seg, dyright_seg, dzright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Acc = [ddxright_seg, ddyright_seg, ddzright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Jerk = [dddxright_seg, dddyright_seg, dddzright_seg ];
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.RBF = [];
            
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.Grp = grpR_seg;
            SegData{gg}.Trial{ii}.Hand{R}.Segment{s}.GrpNum = grp_num;
            
            tskcnt = tskcnt + 1;
        end
        
    end
    All{gg}.taskcount = tskcnt;
    
    %Store all
    PosL = [xleft, yleft, zleft ];
    VelocityL = [dxleft, dyleft, dzleft ];
    AccL = [ddxleft, ddyleft, ddzleft ];
    JerkL = [dddxleft, dddyleft, dddzleft ];
    
    All{gg}.Time = tL;
    All{gg}.Vel = VelocityL;
    All{gg}.Acc = AccL;
    All{gg}.Jerk = JerkL;
    All{gg}.Pos = PosL;
    All{gg}.Grp = grpL ;

end

fprintf('Done \n');

mean(samplength{1})
mean(samplength{2})
mean(samplength{3})

min(samplength{1})
min(samplength{2})
min(samplength{3})

figure
plot(samplength{1})
hold on
plot(samplength{2})
hold on
plot(samplength{3})
hold off


%% Get rbfs for all data DX Vs DDX, DY Vs DDY, DZ Vs DDZ


StrLabels = {'Novice','Expert'};

gamma = 0.6; %0.6
grplist = [1,2,3];

%for printing to CSV for EDGE data viz
filenameviz = {'novice.csv','intermediate.csv','expert.csv'};

fileID = fopen('allgroups.csv','w'); %CSV
fprintf(fileID,'time,id,grp,x,y,z,dx,dy,dz,ddx,ddy,ddz,dddx,dddy,dddz,rbfva1,rbfva2,rbfva3,rbfaj1,rbfaj2,rbfaj3\n');

fprintf('Getting RBFS...');
for gg = grplist
    z1{gg} = [];
    z2{gg} = [];
    z3{gg} = [];
    
%     fileID = fopen(filenameviz{gg},'w'); %CSV
%     fprintf(fileID,'time,id,x,y,z,dx,dy,dz,ddx,ddy,ddz,dddx,dddy,dddz,rbfva1,rbfva2,rbfva3,rbfaj1,rbfaj2,rbfaj3\n');
    gindex = 1;
    
    for ii = 1:length(SegData{gg}.Trial) %limsegs
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand)
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};
                
                
               %Get and sort data for individual grasps (Left)
                if(strcmp(rbfstates,'accjerk'))
                    temp1 = [struct.Acc(:,1), struct.Jerk(:,1)];
                    temp2 = [struct.Acc(:,2), struct.Jerk(:,2)];
                    temp3 = [struct.Acc(:,3), struct.Jerk(:,3)];
                elseif(strcmp(rbfstates,'velacc'))
                    temp1 = [struct.Vel(:,1), struct.Acc(:,1)];
                    temp2 = [struct.Vel(:,2), struct.Acc(:,2)];
                    temp3 = [struct.Vel(:,3), struct.Acc(:,3)];
                elseif(strcmp(rbfstates,'velaccjerk'))
                    temp1 = [struct.Vel(:,1), struct.Acc(:,1), struct.Jerk(:,1)];
                    temp2 = [struct.Vel(:,2), struct.Acc(:,2), struct.Jerk(:,2)];
                    temp3 = [struct.Vel(:,3), struct.Acc(:,3), struct.Jerk(:,3)];
                elseif(strcmp(rbfstates,'posvelacc'))
                    temp1 = [struct.Pos(:,1), struct.Vel(:,1), struct.Acc(:,1)];
                    temp2 = [struct.Pos(:,2), struct.Vel(:,2), struct.Acc(:,2)];
                    temp3 = [struct.Pos(:,3), struct.Vel(:,3), struct.Acc(:,3)];
                else
                    temp1 = [struct.Vel(:,1), struct.Acc(:,1)];
                    temp2 = [struct.Vel(:,2), struct.Acc(:,2)];
                    temp3 = [struct.Vel(:,3), struct.Acc(:,3)];
                end
                rb1 = RadialBasisFunction(temp1,gamma);
                rb2 = RadialBasisFunction(temp2,gamma);
                rb3 = RadialBasisFunction(temp3,gamma);
                
                
                %for additional rbf info
                rt1 = [struct.Acc(:,1), struct.Jerk(:,1)];
                rt2 = [struct.Acc(:,2), struct.Jerk(:,2)];
                rt3 = [struct.Acc(:,3), struct.Jerk(:,3)];
                rbs1 = RadialBasisFunction(rt1,gamma);
                rbs2 = RadialBasisFunction(rt2,gamma);
                rbs3 = RadialBasisFunction(rt3,gamma);
                
                
                %for CSV output
                tout = struct.Time;
                %dat = [tout,ones(length(tout),1)*gindex,struct.Pos,struct.Vel,struct.Acc,struct.Jerk,rb1,rb2,rb3,rbs1,rbs2,rbs3];
                dat = [tout,ones(length(tout),1)*gindex,ones(length(tout),1)*gg,struct.Pos,struct.Vel,struct.Acc,struct.Jerk,rb1,rb2,rb3,rbs1,rbs2,rbs3];
                for ffo = 1:size(dat,1)
                    fprintf(fileID,'%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',dat(ffo,:));
                end
                gindex = gindex + 1;
                
                %stash the RBF back in the struct
                SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.RBF = [rb1, rb2, rb3];

                %for plots
                z1{gg} = [z1{gg}; temp1, rb1];
                z2{gg} = [z2{gg}; temp2, rb2];
                z3{gg} = [z3{gg}; temp3, rb3];

            end
        end
    end
    %fclose(fileID);
end
fprintf('Done \n');

fclose(fileID);


%% getting differences in RBF (cross matrix)

fprintf('Plotting 3D Ratio RBFs...');

gamma = 0.6;
grplist2 = [1,3];

Relative = [];

for gg = grplist2
    Ratio1{gg} = [];
    Ratio2{gg} = [];
    Ratio3{gg} = [];
    for ii = 1:length(SegData{gg}.Trial) %limsegs
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand)
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                
                %get the current group, trial, hand, and segment
                struct1 = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};

                temp1_1 = [struct1.Vel(:,1), struct1.Acc(:,1)];
                temp2_1 = [struct1.Vel(:,2), struct1.Acc(:,2)];
                temp3_1 = [struct1.Vel(:,3), struct1.Acc(:,3)];
                
                %similarities
                D1 = pdist2(temp1_1, temp1_1);
                D2 = pdist2(temp2_1, temp2_1);
                D3 = pdist2(temp3_1, temp3_1);

                D1 = exp(-(D1 .^ 2) ./ ( 2*gamma^2));
                D2 = exp(-(D2 .^ 2) ./ ( 2*gamma^2));
                D3 = exp(-(D3 .^ 2) ./ ( 2*gamma^2));
                
                Sim1 = mean(D1,1)'; %scale by length yo
                Sim2 = mean(D2,1)'; %scale by length yo
                Sim3 = mean(D3,1)'; %scale by length yo

                for gg2 = grplist2
                    Relative{1} = [];
                    Relative{2} = [];
                    Relative{3} = [];
                    if(gg2 ~= gg)
                        for ii2 = 1:length(SegData{gg2}.Trial) %limsegs
                            for hh2 = 1:length(SegData{gg2}.Trial{ii2}.Hand)
                                for ss2 = 1:length(SegData{gg2}.Trial{ii2}.Hand{hh2}.Segment)
                                    %get the current group, trial, hand, and segment
                                    struct2 = SegData{gg2}.Trial{ii2}.Hand{hh2}.Segment{ss2};

                                    temp1_2 = [struct2.Vel(:,1), struct2.Acc(:,1)];
                                    temp2_2 = [struct2.Vel(:,2), struct2.Acc(:,2)];
                                    temp3_2 = [struct2.Vel(:,3), struct2.Acc(:,3)];
                                    
                                    % _1 will be rows, _2 columns
                                    D1 = pdist2(temp1_1, temp1_2);
                                    D2 = pdist2(temp2_1, temp2_2);
                                    D3 = pdist2(temp3_1, temp3_2);

                                    % will be low
                                    D1 = exp(-(D1 .^ 2) ./ ( 2*gamma^2));
                                    D2 = exp(-(D2 .^ 2) ./ ( 2*gamma^2));
                                    D3 = exp(-(D3 .^ 2) ./ ( 2*gamma^2));
                                    
                                    Diff1 = mean(D1,2); %scale by length yo
                                    Diff2 = mean(D2,2); %scale by length yo
                                    Diff3 = mean(D3,2); %scale by length yo

                                    Relative{1} = [Relative{1}, Diff1];
                                    Relative{2} = [Relative{2}, Diff2];
                                    Relative{3} = [Relative{3}, Diff3];
                                end
                            end
                        end
                        
                        Ratio1{gg} = [Ratio1{gg}; temp1_1, Sim1 ./(Sim1 + mean(Relative{1},2) ) ];
                        Ratio2{gg} = [Ratio2{gg}; temp2_1, Sim2 ./(Sim2 + mean(Relative{2},2) ) ];
                        Ratio3{gg} = [Ratio3{gg}; temp3_1, Sim3 ./(Sim3 + mean(Relative{3},2) ) ];
                    end
                end
                
                
                     

            end
        end
    end
end


trAlpha = 0.2;
trPtSize = 1;
trCubes = 1;

figure
gg = 1;
d = [Ratio3{gg}(:,1),Ratio3{gg}(:,2),Ratio3{gg}(:,3)];
h = plot3transparent(d, 'r', trAlpha, trPtSize, trCubes);
hold on
gg = 3;
d = [Ratio3{gg}(:,1),Ratio3{gg}(:,2),Ratio3{gg}(:,3)];
h = plot3transparent(d, 'b', trAlpha*2, trPtSize, trCubes);
hold off
zlim([0 1])
title('Discriminant Pdist')
legend('novices','experts')
xlabel('dx')
ylabel('ddx')
zlabel('rbf1/(rbf1 + rbf2)')

fprintf('Done \n');

%% Try using Adaboost

% http://www.mathworks.com/help/stats/compactclassificationtree.view.html
% http://www.mathworks.com/help/stats/ensemble-methods.html
% http://www.mathworks.com/help/stats/fitensemble.html

XAD = [];
YAD = [];
for gc = [1,3]
    if(size(z1{gc},2) > 3)
        tmp = [z1{gc}(:,4),z2{gc}(:,4),z3{gc}(:,4)];
        XAD = [XAD ; tmp];
        YAD = [YAD ; ones(length(tmp),1)*gc];
    else
        tmp = [z1{gc}(:,3),z2{gc}(:,3),z3{gc}(:,3)];
        XAD = [XAD ; tmp];
        YAD = [YAD ; ones(length(tmp),1)*gc];
    end
end

ClassTreeEns = fitensemble(XAD,YAD,'AdaBoostM1',100,'Tree');

rsLoss = resubLoss(ClassTreeEns,'Mode','Cumulative');

figure
plot(rsLoss);
xlabel('Number of Learning Cycles');
ylabel('Resubstitution Loss');

predY = predict(ClassTreeEns,XAD);

%% Plot RBFS


trAlpha = 0.1;
trPtSize = 1;
trCubes = 1;

fprintf('Plotting 3D RBFs...');

%plot in 3D
figure
gc = 1;
if(size(z1{gc},2) > 3)
    d = [z1{gc}(:,1),z1{gc}(:,2),z1{gc}(:,4)];
else
    d = [z1{gc}(:,1),z1{gc}(:,2),z1{gc}(:,3)];
end
h = plot3transparent(d, 'r', trAlpha, trPtSize, trCubes);
% scatter3(z1{1}(:,1),z1{1}(:,2),z1{1}(:,3),'r*')
hold on
gc = 3;
if(size(z1{gc},2) > 3)
    d = [z1{gc}(:,1),z1{gc}(:,2),z1{gc}(:,4)];
else
    d = [z1{gc}(:,1),z1{gc}(:,2),z1{gc}(:,3)];
end
h = plot3transparent(d, 'b', trAlpha*2, trPtSize, trCubes);
% scatter3(z1{2}(:,1),z1{2}(:,2),z1{2}(:,3),'b+')
hold off
title('3D radial basis x')
legend('novices','experts')
if(strcmp(rbfstates,'velacc'))
    xlabel('dx')
    ylabel('ddx')
elseif(strcmp(rbfstates,'accjerk'))
    xlabel('ddx')
    ylabel('dddx')
elseif(strcmp(rbfstates,'posvelacc'))
    xlabel('x')
    ylabel('dx')
end
zlabel('rbf')

%plot in 3D
figure
gc = 1;
if(size(z2{gc},2) > 3)
    d = [z2{gc}(:,1),z2{gc}(:,2),z2{gc}(:,4)];
else
    d = [z2{gc}(:,1),z2{gc}(:,2),z2{gc}(:,3)];
end
h = plot3transparent(d, 'r', trAlpha, trPtSize, trCubes);
% scatter3(z2{1}(:,1),z2{1}(:,2),z2{1}(:,3),'r*')
hold on
gc = 3;
if(size(z2{gc},2) > 3)
    d = [z2{gc}(:,1),z2{gc}(:,2),z2{gc}(:,4)];
else
    d = [z2{gc}(:,1),z2{gc}(:,2),z2{gc}(:,3)];
end
h = plot3transparent(d, 'b', trAlpha*2, trPtSize, trCubes);
% scatter3(z2{2}(:,1),z2{2}(:,2),z2{2}(:,3),'b+')
hold off
title('3D radial basis y')
legend('novices','experts')
if(strcmp(rbfstates,'velacc'))
    xlabel('dy')
    ylabel('ddy')
elseif(strcmp(rbfstates,'accjerk'))
    xlabel('ddy')
    ylabel('dddy')
elseif(strcmp(rbfstates,'posvelacc'))
    xlabel('y')
    ylabel('dy')
end
zlabel('rbf')

%plot in 3D
figure
gc = 1;
if(size(z3{gc},2) > 3)
    d = [z3{gc}(:,1),z3{gc}(:,2),z3{gc}(:,4)];
else
    d = [z3{gc}(:,1),z3{gc}(:,2),z3{gc}(:,3)];
end
h = plot3transparent(d, 'r', trAlpha, trPtSize, trCubes);
% scatter3(z3{1}(:,1),z3{1}(:,2),z3{1}(:,3),'r*')
hold on
gc = 3;
if(size(z3{gc},2) > 3)
    d = [z3{gc}(:,1),z3{gc}(:,2),z3{gc}(:,4)];
else
    d = [z3{gc}(:,1),z3{gc}(:,2),z3{gc}(:,3)];
end
h = plot3transparent(d, 'b', trAlpha*2, trPtSize, trCubes);
% scatter3(z3{2}(:,1),z3{2}(:,2),z3{2}(:,3),'b+')
hold off
title('3D radial basis z')
legend('novices','experts')
if(strcmp(rbfstates,'velacc'))
    xlabel('dz')
    ylabel('ddz')
elseif(strcmp(rbfstates,'accjerk'))
    xlabel('ddz')
    ylabel('dddz')
elseif(strcmp(rbfstates,'posvelacc'))
    xlabel('z')
    ylabel('dz')
end
zlabel('rbf')

fprintf('Done \n');

%% Plot RBF heights in 3D with bounding box/ hull


trAlpha = 0.2;
trPtSize = 0.01;
trCubes = 1;

figure
gc = 1;
%Get RBF heights and plot those from each dimension
if(size(z1{gc},2) > 3)
    d = [z1{gc}(:,4),z2{gc}(:,4),z3{gc}(:,4)];
else
    d = [z1{gc}(:,3),z2{gc}(:,3),z3{gc}(:,3)];
end
h = plot3transparent(d, 'r', trAlpha, trPtSize, trCubes);
% scatter3(z1{1}(:,3),z2{1}(:,3),z3{1}(:,3),'r*')
hold on
gc = 3;
if(size(z1{gc},2) > 3)
    d = [z1{gc}(:,4),z2{gc}(:,4),z3{gc}(:,4)];
else
    d = [z1{gc}(:,3),z2{gc}(:,3),z3{gc}(:,3)];
end
h = plot3transparent(d, 'b', trAlpha, trPtSize, trCubes);
% scatter3(z1{2}(:,3),z2{2}(:,3),z3{2}(:,3),'b+')
hold on
if(size(z1{gc},2) > 3)
    XT = [z1{gc}(:,4),z2{gc}(:,4),z3{gc}(:,4)];
else
    XT = [z1{gc}(:,3),z2{gc}(:,3),z3{gc}(:,3)];
end
if(strcmp(boundtype,'cube'))
    BoxModel = PlotBoundingCube(XT,0.95,[0 1 0])
elseif(strcmp(boundtype,'hull'))
   [DT,hull] = convHull98Percent(XT,0.98);
   h = convHull3DPlot(DT,hull); 
end
hold on
% PlotQDA3D(MdlQuadratic,D,25);
hold off
title('3D radial basis heights')
legend('novices','experts')
if(strcmp(rbfstates,'velacc'))
    xlabel('rbf(dx,ddx)')
    ylabel('rbf(dy,ddy)')
    zlabel('rbf(dz,ddz)')
elseif(strcmp(rbfstates,'accjerk'))
    xlabel('rbf(ddx,dddx)')
    ylabel('rbf(ddy,dddy)')
    zlabel('rbf(ddz,dddz)')
end

%% Try predicting with bounding cube/ hull for all data (BAD)

mapz = [0,-1,1];
classinit = [];
for ii = grplist
    if(size(z1{ii},2) > 3)
        XT = [z1{ii}(:,4),z2{ii}(:,4),z3{ii}(:,4)];
    else
        XT = [z1{ii}(:,3),z2{ii}(:,3),z3{ii}(:,3)];
    end
    if(strcmp(boundtype,'cube'))
       mask = BoundingBoxCheck(XT,BoxModel);
    elseif(strcmp(boundtype,'hull'))
       mask = checkInsideHull(DT,XT);
    else
       mask = BoundingBoxCheck(XT,BoxModel);
    end
    
    classinit = [classinit; mask, repmat(mapz(ii),length(XT),1) ];
end

acc = classinit(:,1) == classinit(:,2);

sum(acc)/length(acc)


%% classify in cube average for a grasp
hh = 2;

for gg = grplist
    pcb{gg} = [];
    pcgrips{gg} = [];
    for ii = 1:length(SegData{gg}.Trial) % limsegs %
       tempc = [];
       pcgrips{gg}.total{ii} = [];
       for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
           %Get RBF for individual grasps
            XT = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.RBF;
            
            if(strcmp(boundtype,'cube'))
               classification = BoundingBoxCheck(XT,BoxModel);
            elseif(strcmp(boundtype,'hull'))
               classification = checkInsideHull(DT,XT);
            else
               classification = BoundingBoxCheck(XT,BoxModel);
            end

            classification = ~classification;
            amtout = sum(classification)/length(classification);
            
            pcgrips{gg}.total{ii}(ss) = amtout;
            

            tempc = [tempc; amtout];
       end
       
       pcb{gg}(ii) = mean(tempc);
       
       figure
       plot(pcgrips{gg}.total{ii})
       str = sprintf('amount outside of box for class %i, user %i', gg, ii);
       title(str);
       xlabel('segment #')
       ylabel('total outside box')
       ylim([0 1]) 
    end
end

figure
scatter(ones(length(pcb{1}),1), pcb{1},'r*')
hold on
scatter(ones(length(pcb{3}),1)*2, pcb{3},'b+')
hold off
title('average outside of box')
legend('novice','expert')
xlabel('segment #')
ylabel('total outside box')

figure
plot(pcb{1})
title('total outside per user novices')
xlabel('user #')
ylabel('total outside box')
axis([0 30 0 1])

figure
plot(pcb{3})
title('total outside per user experts')
xlabel('user #')
ylabel('total outside box')
axis([0 6 0 1])

%% LEAVE ONE OUT CLASSIFY

grplist = [1,2,3];

gnov = 1;
gint = 2;
gexp = 3;

SepHands = true;

BoxModelLOO = [];
DTLOO = [];
hullLOO = [];

storeIntmean = [];
storeEXPAll = [];
storeEXPmean = [];
storeEXPstd = [];
storeNOVmean = [];

fprintf('Performing leave one out...');

for e_o = 1:length(SegData{gexp}.Trial)  %loop through all expert surgeons (leave out)
    e_o
    %train box using leave one out experts
    XGE = getRBFdata(SegData,gexp,e_o,SepHands);

    %compute bounding limits
    if(SepHands)
        if(strcmp(boundtype,'cube'))
            BoxModelLOO{1} = PlotBoundingCube(XGE{1},1,[0 1 0],0);
            BoxModelLOO{2} = PlotBoundingCube(XGE{2},1,[0 1 0],0);
        elseif(strcmp(boundtype,'hull'))
           [DTLOO{1},hullLOO{1}] = convHull98Percent(XGE{1},0.98);
           [DTLOO{2},hullLOO{2}] = convHull98Percent(XGE{2},0.98);
        end
    else
        if(strcmp(boundtype,'cube'))
            BoxModelLOO = PlotBoundingCube(XGE,1,[0 1 0],0);
        elseif(strcmp(boundtype,'hull'))
           [DTLOO,hullLOO] = convHull98Percent(XGE,0.98);
        end
    end
    
    
    %online for experts
   for ii = 1:length(SegData{gexp}.Trial) %loop through all novice surgeons keep some
       stashAllExp = [];
       for hh = 1:length(SegData{gexp}.Trial{ii}.Hand)
           for ss = 1:length(SegData{gexp}.Trial{ii}.Hand{hh}.Segment) %all segments from that surgeon

                %put rbfs in big old matrix
                XGEO =  SegData{gexp}.Trial{ii}.Hand{hh}.Segment{ss}.RBF;
    
                %check classification
                if(SepHands)
                    if(strcmp(boundtype,'cube'))
                       class = BoundingBoxCheck(XGEO,BoxModelLOO{hh});
                    elseif(strcmp(boundtype,'hull'))
                       class = checkInsideHull(DTLOO{hh},XGEO);
                    else
                       class = BoundingBoxCheck(XGEO,BoxModelLOO{hh});
                    end
                else
                    if(strcmp(boundtype,'cube'))
                       class = BoundingBoxCheck(XGEO,BoxModelLOO);
                    elseif(strcmp(boundtype,'hull'))
                       class = checkInsideHull(DTLOO,XGEO);
                    else
                       class = BoundingBoxCheck(XGEO,BoxModelLOO);
                    end
                end

                class = ~class;
                
                nrm = NormRowWise(XGEO);
                
                stashAllExp = [ stashAllExp; class];

           end
       end
       if(ii == e_o)
            storeEXPAll{e_o} = stashAllExp;
            storeEXPmean(e_o) = mean(stashAllExp);
            storeEXPstd(e_o) = std(stashAllExp);
       end
   end
   
   %get intermediates with leave one out
   for ii = 1:length(SegData{gint}.Trial) %loop through all novice surgeons keep some
       %only use surgeon who is not being left out
       stashAllInt = [];
       for hh = 1:length(SegData{gint}.Trial{ii}.Hand)
           for ss = 1:length(SegData{gint}.Trial{ii}.Hand{hh}.Segment) %all segments from that surgeon

                %put rbfs in big old matrix
                XGI =  SegData{gint}.Trial{ii}.Hand{hh}.Segment{ss}.RBF;

                if(SepHands)
                    if(strcmp(boundtype,'cube'))
                       classification = BoundingBoxCheck(XGI,BoxModelLOO{hh});
                    elseif(strcmp(boundtype,'hull'))
                       classification = checkInsideHull(DTLOO{hh},XGI);
                    else
                       classification = BoundingBoxCheck(XGI,BoxModelLOO{hh});
                    end
                else
                    if(strcmp(boundtype,'cube'))
                       classification = BoundingBoxCheck(XGI,BoxModelLOO);
                    elseif(strcmp(boundtype,'hull'))
                       classification = checkInsideHull(DTLOO,XGI);
                    else
                       classification = BoundingBoxCheck(XGI,BoxModelLOO);
                    end
                end

                classification = ~classification;
                nrm = NormRowWise(XGI);
                
                stashAllInt = [ stashAllInt; classification];
           end
       end
       storeIntmean(e_o,ii) = mean(stashAllInt);
   end
   
   %get novices with leave one out
   for ii = 1:length(SegData{gnov}.Trial) %loop through all novice surgeons keep some
       %only use surgeon who is not being left out
       stashAllNov = [];
       for hh = 1:length(SegData{gnov}.Trial{ii}.Hand)
           for ss = 1:length(SegData{gnov}.Trial{ii}.Hand{hh}.Segment) %all segments from that surgeon

                %put rbfs in big old matrix
                XGN =  SegData{gnov}.Trial{ii}.Hand{hh}.Segment{ss}.RBF;

                if(SepHands)
                    if(strcmp(boundtype,'cube'))
                       classification = BoundingBoxCheck(XGN,BoxModelLOO{hh});
                    elseif(strcmp(boundtype,'hull'))
                       classification = checkInsideHull(DTLOO{hh},XGN);
                    else
                       classification = BoundingBoxCheck(XGN,BoxModelLOO{hh});
                    end
                else
                    if(strcmp(boundtype,'cube'))
                       classification = BoundingBoxCheck(XGN,BoxModelLOO);
                    elseif(strcmp(boundtype,'hull'))
                       classification = checkInsideHull(DTLOO,XGN);
                    else
                       classification = BoundingBoxCheck(XGN,BoxModelLOO);
                    end
                end

                classification = ~classification;
                nrm = NormRowWise(XGN);
                
                stashAllNov = [ stashAllNov; classification];
           end
       end
       storeNOVmean(e_o,ii) = mean(stashAllNov);
   end
   
   
end

fprintf('Done \n');

[valnov,inov] = sort(mean(storeNOVmean));
[valexp,iexp] = sort(storeEXPmean,'descend');
storeNOVstd = std(storeNOVmean);

figure
errorbar(1:length(valnov),100*valnov,100*storeNOVstd(inov))
hold on
% figure
% plot(1:length(valexp),100*valexp)
errorbar(1:length(valexp),100*valexp,100*storeEXPstd(iexp))
hold off
str = sprintf('Grasp outside box per user (Task: %s)',TaskStr{tsk});
title(str)
xlabel('Surgeon #')
ylabel('% of grasp outside box')
legend('Novices','Experts')


%% plot box score vs fls score

%get the leave one out scores for each trial
boxscores{1} = mean(storeNOVmean);
boxscores{2} = mean(storeIntmean);
boxscores{3} = storeEXPmean;

ignoreint = false;

%pack the fls and box score into one matrix
Correlate = [];
boxall = [];
ggbox = [];
for gg = grplist
    
    if(ignoreint)
        if(gg ~= 2)
            for ii = 1:length(SegData{gg}.Trial) % limsegs
                Correlate{gg}(ii,:) = [SegData{gg}.Trial{ii}.FLSScore, boxscores{gg}(ii) ];
            end

            %get a fit of the line
            boxall = [boxall; boxscores{gg}' ]; %
            ggbox = [ggbox; ones(length( boxscores{gg}),1)*gg ];

            %Correlate{gg}(Correlate{gg}(:,1) <= 0,:) = [];
        end
    else
        for ii = 1:length(SegData{gg}.Trial) % limsegs
            Correlate{gg}(ii,:) = [SegData{gg}.Trial{ii}.FLSScore, boxscores{gg}(ii) ];
        end

        %get a fit of the line
        boxall = [boxall; boxscores{gg}' ]; %
        ggbox = [ggbox; ones(length( boxscores{gg}),1)*gg ];

        %Correlate{gg}(Correlate{gg}(:,1) <= 0,:) = [];
    end
end

%cost matrix (to make it harder to mis classify experts)
costmat = [0, 1, 2; %novices
           1, 0, 2; %intermediates
           5, 5, 0]; %experts

%train a linear classifier
lda = fitcdiscr(boxall,ggbox,'Cost',costmat);

%try predicting
[predY,score] = predict(lda,boxall);
corr = ggbox == predY;
disp('accuracy: ')
acclda = sum(corr)/length(corr)

%borders from lda
brd12 = -lda.Coeffs(1, 2).Const / lda.Coeffs(1, 2).Linear;
brd23 = -lda.Coeffs(2, 3).Const / lda.Coeffs(2, 3).Linear;

%get a fit of the line
if(ignoreint)
    correlateall = [Correlate{1}; Correlate{3}]; %
else
    correlateall = [Correlate{1}; Correlate{2}; Correlate{3}]; %
end

xc1 = correlateall(:,1);
xc2 = correlateall(:,2);
[p,S] = polyfit(xc1,xc2,3)
yCalc = polyval(p,sort(xc1));
[r2 rmse] = rsquare(xc2,yCalc)

%plot fls vs dockter box score
figure
gc = 1;
scatter(Correlate{gc}(:,1),Correlate{gc}(:,2),'rx');
hold on
gc = 2;
if(~ignoreint)
    scatter(Correlate{gc}(:,1),Correlate{gc}(:,2),'cd');
    hold on
end
gc = 3;
scatter(Correlate{gc}(:,1),Correlate{gc}(:,2),'b*');
hold on
plot(sort(xc1),yCalc,'g--') %the fit
hold on
plot([min(correlateall(:,1)), max(correlateall(:,1))], [brd12, brd12], 'c-'); %boundaries
hold on
plot([min(correlateall(:,1)), max(correlateall(:,1))], [brd23, brd23], 'b-'); %boundaries
hold off
str = sprintf('Box Score vs. FLS (Task: %s, RBF: %s, Gamma: %f)',TaskStr{tsk}, rbfstates, gamma);
title(str)
xlabel('fls score')
ylabel('dockter-box score')
legend('novices','Intermed','experts','fit','LDA 1-2','LDA 2-3')


%Plot classification via lda
figure
gscatter(correlateall(:,1),correlateall(:,2),predY)
xlabel('fls score')
ylabel('dockter-box score')
legend('novices','Intermed','experts')
strt = sprintf('classification via lda, accuracy = %f', acclda);
title(strt)



%% Now try actually gradings

%Use demerit system where total time spent outside box for a given segment
%counts against your score

grplist = [1,3];

for gg = grplist
    userGrade{gg} = [];
    for ii = 1:length(SegData{gg}.Trial) % limsegs %
       demerit = 0;
       for hh = 1:length(SegData{gg}.Trial{ii}.Hand)
           for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)

                XT = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.RBF;

                if(strcmp(boundtype,'cube'))
                   classification = BoundingBoxCheck(XT,BoxModel);
                elseif(strcmp(boundtype,'hull'))
                   classification = checkInsideHull(DT,XT);
                else
                   classification = BoundingBoxCheck(XT,BoxModel);
                end
                classification = ~classification;
                amtout = sum(classification)/length(classification);

                demerit = demerit + amtout*100;

           end
       
       end
       startgrade = 100 - ((2*demerit)/(ss));
       
       userGrade{gg}(ii) = startgrade;
    end
end

figure
plot(userGrade{1},'r-')
hold on
plot(userGrade{2},'g-')
hold on
plot(userGrade{3},'b-')
hold off
title('final grade')
axis([0 29 0 100])


return


%% Now Get RBFS limiting time segments (DONT USE)

limsegs = min([length(SegData{1}.Trial), length(SegData{3}.Trial)]);

grplist = [1,3];

StrLabels = {'Novice','Expert'};
labelsall = [];

gamma = 2; %RBF func

tlim = 1; %use this to set wndow size in data

hh = 1;
for gg = grplist
    z1s{gg} = [];
    z2s{gg} = [];
    z3s{gg} = [];
    for ii = 1:length(SegData{gg}.Trial) %limsegs
       for ss = 1:length(SegData{gg}.Trial{ii}.Hand{1}.Segment)
           segsize = length(SegData{gg}.Trial{ii}.Hand{1}.Segment{ss}.Time);
           if(segsize > tlim)
               %only grab the end of the grasp
               idcs = 1:segsize; 
               %idcs = segsize-tlim:segsize; %ending
               %idcs = 1:tlim; %beggining
               
               %Get and sort data for individual grasps
                if(strcmp(rbfstates,'accjerk'))
                    temp1 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,1), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Jerk(idcs,1)];
                    temp2 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,2), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Jerk(idcs,2)];
                    temp3 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,3), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Jerk(idcs,3)];
                elseif(strcmp(rbfstates,'velacc'))
                    temp1 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Vel(idcs,1), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,1)];
                    temp2 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Vel(idcs,2), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,2)];
                    temp3 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Vel(idcs,3), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,3)];
                else
                    temp1 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Vel(idcs,1), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,1)];
                    temp2 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Vel(idcs,2), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,2)];
                    temp3 = [SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Vel(idcs,3), SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.Acc(idcs,3)];
                end
                z1temp = RadialBasisFunction(temp1,gamma);
                z2temp = RadialBasisFunction(temp2,gamma);
                z3temp = RadialBasisFunction(temp3,gamma);
                
                SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.RBF = [z1temp, z2temp, z3temp];

                z1s{gg} = [z1s{gg}; z1temp ];
                z2s{gg} = [z2s{gg}; z2temp ];
                z3s{gg} = [z3s{gg}; z3temp ];

                labelsall = [ labelsall; SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.GrpNum];
           end
       end
    end
end


trAlpha = 0.3;
trPtSize = 0.01;
trCubes = 1;

%plot in 3D
figure
gc = 1;
d = [z1s{gc},z2s{gc},z3s{gc}];
h = plot3transparent(d, 'r', trAlpha, trPtSize, trCubes);
% scatter3(z1s{1},z2s{1},z3s{1},'r*')
hold on
gc = 3;
d = [z1s{gc},z2s{gc},z3s{gc}];
h = plot3transparent(d, 'b', trAlpha, trPtSize*2, trCubes);
% scatter3(z1s{3},z2s{3},z3s{3},'b+')
hold on
XT = [z1s{3},z2s{3},z3s{3}];
if(strcmp(boundtype,'cube'))
    BoxModel = PlotBoundingCube(XT,1,[0 1 0]) %Get the cube!!!!!!!!!!!!!!!!!!!!!!
elseif(strcmp(boundtype,'hull'))
   [DT,hull] = convHull98Percent(XT,0.98);
   h = convHull3DPlot(DT,hull); 
end

hold off
title('3D radial basis heights')
legend('novices','experts')
if(strcmp(rbfstates,'velacc'))
    xlabel('rbf(dx,ddx)')
    ylabel('rbf(dy,ddy)')
    zlabel('rbf(dz,ddz)')
elseif(strcmp(rbfstates,'accjerk'))
    xlabel('rbf(ddx,dddx)')
    ylabel('rbf(ddy,dddy)')
    zlabel('rbf(ddz,dddz)')
end

XTc = [z1s{1},z2s{1},z3s{1}];
if(strcmp(boundtype,'cube'))
   mask = BoundingBoxCheck(XTc,BoxModel);
elseif(strcmp(boundtype,'hull'))
   mask = checkInsideHull(DT,XTc);
else
   mask = BoundingBoxCheck(XTc,BoxModel);
end

sum(mask)/length(mask)

%% OLD Train matlabs svm

counts = [];
jj=1;%nov
counts = [counts, length(All{jj}.Acc) ];
jj=3;%exp
counts = [counts, length(All{jj}.Acc) ];

maxes  = min( counts);
limz= 1:maxes;

X=[];
jj=1;%nov
X = [X;  All{jj}.Vel(limz,1), All{jj}.Vel(limz,2), All{jj}.Vel(limz,3), All{jj}.Acc(limz,1), All{jj}.Acc(limz,2), All{jj}.Acc(limz,3)  ];
jj=3;%exp
X = [X;  All{jj}.Vel(limz,1), All{jj}.Vel(limz,2), All{jj}.Vel(limz,3), All{jj}.Acc(limz,1), All{jj}.Acc(limz,2), All{jj}.Acc(limz,3)  ];


Y = [ones(maxes,1)*1; ones(maxes,1)*3 ]; %nov = 0, exp = 1


SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');

[est,score] = predict(SVMModel,X);

plot(est)

corr = est == Y;
acc = sum(corr)/length(corr)

%% test each matlab svm segment



grplist = [1,3];

for gg = grplist
    id = 1;
    
    for ii = 1:length(SegData{gg}.Trial) 
       svmTrial{gg}.trial{ii} = []; 
       
       for ss = 1:length(SegData{gg}.Trial{ii}.Hand{1}.Segment)
           %Get and sort data for individual grasps
           goop = SegData{gg}.Trial{ii}.Hand{1}.Segment{ss};
           tempx = [goop.Vel(:,1), goop.Vel(:,2), goop.Vel(:,3), goop.Acc(:,1), goop.Acc(:,2), goop.Acc(:,3)];
       
           [est,score] = predict(SVMModel,tempx);
           
           svmTrial{gg}.trial{ii} = [svmTrial{gg}.trial{ii}; score];
            
           svmsegmentEst{gg}(id,:) = mode(est);
           svmsegmentScore{gg}(id,:) = sum(score,1);
           id = id +1;
            
       end
       
       svmTrial{gg}.totalscore{ii} = sum(svmTrial{gg}.trial{ii},1); 
       [~,id] = max(svmTrial{gg}.totalscore{ii});
       svmTrial{gg}.class{ii} = grplist(id);
       
       figure
       plot(svmTrial{gg}.trial{ii})
       str = sprintf('class %i, trial %i',gg, ii);
       title(str);
       legend('nov','exp');
    end
end

corr1 = svmsegmentEst{1} == 1;
corr3 = svmsegmentEst{3} == 3;

acc1 = sum(corr1)/length(corr1)
acc3 = sum(corr3)/length(corr3)


return