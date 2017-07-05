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

%% first load in this large data structure, if isnt already

if( ~exist('DataGlb'))
    fprintf('Loading Global Data Structure EdgeDataGlb.mat (huge) ...');
    try
        load('D:\temp\EdgeDataGlb.mat')
    catch
        disp('Could not load local. Loading from M drive...');
        load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EdgeDataGlb.mat')
    end
    try
        load('D:\temp\EDGE_Segments.mat')
    catch
        disp('Could not load local. Loading from M drive...');
        load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EDGE_Segments.mat')
    end
    fprintf(' DONE.\n');
end

tsk = 1; % Peg Transfer task
myGroups = [g.gtExp g.flsInt g.flsNov ]; % select only three skill groups

% DEMO CODE for Segments:
% to get the start and end of a given segment: do this:
% Choose Task.
tsk =1 ; %PegTx
% Chose Log.
i = 4;% my specific LogIdx (get this from grp scructure for experts, novices, etc and for loop it
% Choose Segment Type
SegType = SegS.betweenFgActs; % (based only on Forces, excludes position
SegType = SegS.betweenGrActs; % (includes forces and grasper position)

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
          DataGlb.dataLog{i}( iSamples, G.toolPathL),'.b'); hold on;

end


%% Dump segments to text file


GG = [randi(4,10,1),rand(10,1),ones(10,1) ];

fileouts = {'EDGE_segment_novice.csv'; 'EDGE_segment_intermediate.csv'; 'EDGE_segment_expert.csv'};
find = 1;


myGroups = fliplr([g.gtExp g.flsInt g.flsNov ]);


for gg = 1:length(myGroups) %looping through array of logs
    % Initialize
    t = [];
    xleft = [];
    yleft = [];
    zleft = [];
    input_x = [];
    input_y = [];
    input_z = [];
    structData = [];
    outtxt = [];
    
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)})
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        % ADD SEGMENTS:
        % Extract actual segment data: (left first)
        mySegmentsL = SegScheme{tsk}.dataLogSegs(i, SegType ,L); % Left hand
        
        for s = 1:length(mySegmentsL.sgidx) %loop through all segments for that particular grasp   
            
            % selects only the samples of this segment
            iSamples = mySegmentsL.sgidx{s};
            iEnd = iSamples(end);
            iStart = iSamples(1);
            
            t = [ t ; DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iEnd, G.Time)] ;
            xleft = [ xleft ; DataGlb.dataLog{i}( iSamples, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL)] ;
            yleft = [ yleft ; DataGlb.dataLog{i}( iSamples, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL)] ;
            zleft = [ zleft ; DataGlb.dataLog{i}( iSamples, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL)] ;
            

            %take input as final location of motion
            input_x = [input_x; (DataGlb.dataLog{i}( iStart, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL))*ones(length(iSamples),1) ];
            input_y = [input_y; (DataGlb.dataLog{i}( iStart, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL))*ones(length(iSamples),1) ];
            input_z = [input_z; (DataGlb.dataLog{i}( iStart, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL))*ones(length(iSamples),1) ];
            
            
            %This is dumb
            structData = [structData ; DataGlb.dataLog{i}( iSamples,2:end ) ];

        end
        
        
        
    end
    
    
    %output text matrix
    outtxt = [t, xleft, yleft, zleft, structData , input_x, input_y, input_z, gg*ones(length(t),1) ];
    
    %csvwrite(fileouts{find},outtxt)
    
    fid = fopen(fileouts{find}, 'w');
    
    fprintf(fid, 'time,x_corr,y_corr,z_corr,Q1L,Q2L,dL,Q3L,QgL,FgL,Q1R,Q2R,dR,Q3R,QgR,FgR,xL,yL,zL,xR,yR,zR,toolPathL,toolPathR,toolSpeedL,toolSpeedR,toolAccL,toolAccR,toolJerkL,toolJerkR,toolPathCrvL,toolPathCrvR,dQ3L,dQgL,dFgL,dQ3R,dQgR,dFgR,dQ3Latan,dQ3Ratan,dxL,dyL,dzL,dxR,dyR,dzR,inputx,inputy,inputz,skillLevel\n');
    fclose(fid);

    dlmwrite(fileouts{find}, outtxt, '-append', 'precision', '%.6f', 'delimiter', ',');
    find = find + 1;
    
   
end



%% Transform each segment to a new coordinate system: Origin is the end of the segment
ggLbl = {'r', 'g', 'b'};

myGroups = fliplr([g.gtExp g.flsInt g.flsNov ]);

All = []
All{1}.Input = [];
All{1}.Speed = [];
All{1}.Acc = [];
All{1}.Pos  = [];
All{1}.Grp = [];

All{2}.Input = [];
All{2}.Speed = [];
All{2}.Acc = [];
All{2}.Pos  = [];
All{2}.Grp = [];

All{3}.Input = [];
All{3}.Speed = [];
All{3}.Acc = [];
All{3}.Pos  = [];
All{3}.Grp = [];

Online = []
Online{1}.Input = [];
Online{1}.Speed = [];
Online{1}.Acc = [];
Online{1}.Pos  = [];
Online{1}.Grp = [];

Online{2}.Input = [];
Online{2}.Speed = [];
Online{2}.Acc = [];
Online{2}.Pos  = [];
Online{2}.Grp = [];

Online{3}.Input = [];
Online{3}.Speed = [];
Online{3}.Acc = [];
Online{3}.Pos  = [];
Online{3}.Grp = [];

for gg = 1:length(myGroups) %looping through array of logs
    % Initialize
    t = [];
    xleft = [];
    yleft = [];
    zleft = [];
    input_x = [];
    input_y = [];
    input_z = [];
    dxleft = [];
    dyleft = [];
    dzleft = [];
    ddxleft = [];
    ddyleft = [];
    ddzleft = [];
    grpL = [];

    t_o = [];
    xleft_o = [];
    yleft_o = [];
    zleft_o = [];
    input_x_o = [];
    input_y_o = [];
    input_z_o = [];
    dxleft_o = [];
    dyleft_o = [];
    dzleft_o = [];
    ddxleft_o = [];
    ddyleft_o = [];
    ddzleft_o = [];
    grpL_o = [];
    
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)})
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        % ADD SEGMENTS:
        % Extract actual segment data: (left first)
        mySegmentsL = SegScheme{tsk}.dataLogSegs(i, SegType ,L); % Left hand
        
        for s = 1:length(mySegmentsL.sgidx) %loop through all segments for that particular grasp   
            
            % selects only the samples of this segment
            iSamples = mySegmentsL.sgidx{s};
            iEnd = iSamples(end);
            iStart = iSamples(1);
            
            if(ii == length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)}) )
                t_o = [ t_o ; DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iEnd, G.Time)] ;
                xleft_o = [ xleft_o ; DataGlb.dataLog{i}( iSamples, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL)] ;
                yleft_o = [ yleft_o ; DataGlb.dataLog{i}( iSamples, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL)] ;
                zleft_o = [ zleft_o ; DataGlb.dataLog{i}( iSamples, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL)] ;

                %take input as final location of motion
                input_x_o = [input_x_o; (DataGlb.dataLog{i}( iStart, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL))*ones(length(iSamples),1) ];
                input_y_o = [input_y_o; (DataGlb.dataLog{i}( iStart, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL))*ones(length(iSamples),1) ];
                input_z_o = [input_z_o; (DataGlb.dataLog{i}( iStart, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL))*ones(length(iSamples),1) ];

                dxleft_o = [ dxleft_o ; DataGlb.dataLog{i}( iSamples, G.dxL) ] ;
                dyleft_o = [ dyleft_o ; DataGlb.dataLog{i}( iSamples, G.dyL) ] ;
                dzleft_o = [ dzleft_o ; DataGlb.dataLog{i}( iSamples, G.dzL) ] ;

                grpL_o = [ grpL_o; repmat( DataGlb.grp.tag(myGroups(gg)), [length(iSamples),1]) ];
                
            else
            
                t = [ t ; DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iEnd, G.Time)] ;
                xleft = [ xleft ; DataGlb.dataLog{i}( iSamples, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL)] ;
                yleft = [ yleft ; DataGlb.dataLog{i}( iSamples, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL)] ;
                zleft = [ zleft ; DataGlb.dataLog{i}( iSamples, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL)] ;

                %take input as final location of motion
                input_x = [input_x; (DataGlb.dataLog{i}( iStart, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL))*ones(length(iSamples),1) ];
                input_y = [input_y; (DataGlb.dataLog{i}( iStart, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL))*ones(length(iSamples),1) ];
                input_z = [input_z; (DataGlb.dataLog{i}( iStart, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL))*ones(length(iSamples),1) ];

                dxleft = [ dxleft ; DataGlb.dataLog{i}( iSamples, G.dxL) ] ;
                dyleft = [ dyleft ; DataGlb.dataLog{i}( iSamples, G.dyL) ] ;
                dzleft = [ dzleft ; DataGlb.dataLog{i}( iSamples, G.dzL) ] ;

                grpL = [ grpL; repmat( DataGlb.grp.tag(myGroups(gg)), [length(iSamples),1]) ];
            end
       
        end
        
        if(ii == length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)}) )
            
            T_o = mean(diff(t_o));
            ddxleft_o = Calculate_velocity( dxleft_o, T_o, 'holobrodko') ;
            ddyleft_o = Calculate_velocity( dyleft_o, T_o, 'holobrodko') ;
            ddzleft_o = Calculate_velocity( dzleft_o, T_o, 'holobrodko') ;


            InputL_o = [input_x_o, input_y_o, input_z_o ];
            PosL_o = [xleft_o, yleft_o, zleft_o ];
            SpeedL_o = [dxleft_o, dyleft_o, dzleft_o ];
            AccL_o = [ddxleft_o, ddyleft_o, ddzleft_o ];
    
            Online{gg}.Input = InputL_o;
            Online{gg}.Speed = SpeedL_o;
            Online{gg}.Acc = AccL_o;
            Online{gg}.Pos  =  PosL_o;
            Online{gg}.Grp   = grpL_o;
                
        end
    end
    
    T = mean(diff(t));
    ddxleft = Calculate_velocity( dxleft, T, 'holobrodko') ;
    ddyleft = Calculate_velocity( dyleft, T, 'holobrodko') ;
    ddzleft = Calculate_velocity( dzleft, T, 'holobrodko') ;
    
    
    InputL = [input_x, input_y, input_z ];
    PosL = [xleft, yleft, zleft ];
    SpeedL = [dxleft, dyleft, dzleft ];
    AccL = [ddxleft, ddyleft, ddzleft ];
    

    All{gg}.Input = InputL;
    All{gg}.Speed = SpeedL;
    All{gg}.Acc = AccL;
    All{gg}.Pos  =  PosL;
    All{gg}.Grp   = grpL ;
    

end


%% Rods KDEs

jj=1;%nov
dat_nov = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=2;%int
dat_int = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=3;%exp
dat_exp = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];

jj=1;%nov
dat_nov_online = [abs( Online{jj}.Pos(:,1)+Online{jj}.Pos(:,2)+Online{jj}.Pos(:,3) ) , abs( Online{jj}.Speed(:,1)+Online{jj}.Speed(:,2)+Online{jj}.Speed(:,3) ), abs( Online{jj}.Acc(:,1)+Online{jj}.Acc(:,2)+Online{jj}.Acc(:,3) ) ];
jj=2;%int
dat_int_online = [abs( Online{jj}.Pos(:,1)+Online{jj}.Pos(:,2)+Online{jj}.Pos(:,3) ) , abs( Online{jj}.Speed(:,1)+Online{jj}.Speed(:,2)+Online{jj}.Speed(:,3) ), abs( Online{jj}.Acc(:,1)+Online{jj}.Acc(:,2)+Online{jj}.Acc(:,3) ) ];
jj=3;%exp
dat_exp_online = [abs( Online{jj}.Pos(:,1)+Online{jj}.Pos(:,2)+Online{jj}.Pos(:,3) ) , abs( Online{jj}.Speed(:,1)+Online{jj}.Speed(:,2)+Online{jj}.Speed(:,3) ), abs( Online{jj}.Acc(:,1)+Online{jj}.Acc(:,2)+Online{jj}.Acc(:,3) ) ];


dat = [dat_nov; dat_int; dat_exp];

figure
scatter3(dat_nov(:,1),dat_nov(:,2),dat_nov(:,3),10,'r.')
hold on
scatter3(dat_int(:,1),dat_int(:,2),dat_int(:,3),10,'b.')
hold on
scatter3(dat_exp(:,1),dat_exp(:,2),dat_exp(:,3),10,'g.')
hold off

title('Training Distributions')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Novice','Intermediate','Expert')

%%
% compute distro
bw = 0.8;
%pdf_nov = kde3d(dat_nov,dat,'graph');
%pdf_int = kde3d(dat_int,dat,'graph');
%tic
%pdf_exp = kde3d(dat_exp,dat,'full');
%toc

%% online (TIM IDEA 7_23_2015)
dat_nov_s = subSampleVec(dat_nov,0.25);
dat_int_s = subSampleVec(dat_int,0.5);
dat_exp_s = subSampleVec(dat_exp,1);

samp = dat_nov_online; %datasample(dat_nov_s,1);

on_class = [0,0,0];

for kk = 1:length(samp)
    disp(sprintf('Loop %i of %i \n', kk, length(samp) ));
    onk = samp(kk,:);

    f_n = kde3d_online(dat_nov_s,onk);
    f_i = kde3d_online(dat_int_s,onk);
    f_e = kde3d_online(dat_exp_s,onk);
    
    est(kk,:) = [f_n,f_i,f_e];
    %[m,idx] = max(est);
    %on_class(idx) = on_class(idx)  + 1;
    on_class = on_class + est(kk,:);
end

[max,idx] = max(on_class)

%% Scale

% pdf_nov.prob = pdf_nov.prob/max(pdf_nov.prob);
% pdf_int.prob = pdf_int.prob/max(pdf_int.prob);
% pdf_exp.prob = pdf_exp.prob/max(pdf_exp.prob);

%%
figure
scatter3(dat_nov(:,1),dat_nov(:,2),dat_nov(:,3),8,pdf_nov.prob);
colormap(cool);
colorbar;
title('Density Estimates Novice')
xlabel('pos')
ylabel('velocity')
zlabel('acceleration')
% 
% figure
% scatter3(dat_int(:,1),dat_int(:,2),dat_int(:,3),8,pdf_int.prob);
% colormap(cool);
% colorbar;
% title('Density Estimates Intermediate')
% xlabel('pos')
% ylabel('velocity')
% zlabel('acceleration')
% 
% 
% figure
% scatter3(dat_exp(:,1),dat_exp(:,2),dat_exp(:,3),8,pdf_exp.prob);
% colormap(cool);
% colorbar;
% title('Density Estimates Expert')
% xlabel('pos')
% ylabel('velocity')
% zlabel('acceleration')
% %%
% [class_nov,sums_nov] = SkillClassificationPDF(dat_nov(:,[1,2,3]),pdf_nov,pdf_int,pdf_exp)
% 
% [class_int,sums_int] = SkillClassificationPDF(dat_int(:,[1,2,3]),pdf_nov,pdf_int,pdf_exp)
% 
% [class_exp,sums_exp] = SkillClassificationPDF(dat_exp(:,[1,2,3]),pdf_nov,pdf_int,pdf_exp)
%                                                                                                                                                                                                                                                     
return