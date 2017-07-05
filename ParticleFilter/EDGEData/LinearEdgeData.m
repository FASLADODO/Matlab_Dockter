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
            
            dxleft = [ dxleft ; DataGlb.dataLog{i}( iSamples, G.dxL) - DataGlb.dataLog{i}( iEnd, G.dxL)] ;
            dyleft = [ dyleft ; DataGlb.dataLog{i}( iSamples, G.dyL) - DataGlb.dataLog{i}( iEnd, G.dyL)] ;
            dzleft = [ dzleft ; DataGlb.dataLog{i}( iSamples, G.dzL) - DataGlb.dataLog{i}( iEnd, G.dzL)] ;
            
            grpL = [ grpL; repmat( DataGlb.grp.tag(myGroups(gg)), [length(iSamples),1]) ];
       
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
    
    if(ii == length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)}) )
        Online{gg}.Input = InputL;
        Online{gg}.Speed = SpeedL;
        Online{gg}.Acc = AccL;
        Online{gg}.Pos  =  PosL;
        Online{gg}.Grp   = grpL ;
    else
        All{gg}.Input = InputL;
        All{gg}.Speed = SpeedL;
        All{gg}.Acc = AccL;
        All{gg}.Pos  =  PosL;
        All{gg}.Grp   = grpL ;
    end
    

end

%% Plot every dimension for fun

dirstr = {'x', 'y', 'z'};
dirnum = [1,2,3];

kk = 1;
for ii = 1:3 %x,y,z
    for  jj = 1:length(myGroups)
        figure(kk)

        %dist = abs( All{jj}.Pos(:,1) + All{jj}.Pos(:,2) + All{jj}.Pos(:,3) );
        %vel = abs( All{jj}.Speed(:,1) + All{jj}.Speed(:,2) + All{jj}.Speed(:,3) );

        scatter3(All{jj}.Pos(:,ii),All{jj}.Speed(:,ii),All{jj}.Acc(:,ii),'r');

        strtitle = strcat(DataGlb.grp.tag(myGroups(jj)),' phase portrait in: ',dirstr(ii));
        title(strtitle)
        xlabel('position')
        ylabel('velocity')
        zlabel('acceleration')
        
        kk = kk + 1;
    end
end

%% Gaussian pdfs

jj=1;%nov
X1 = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=2;%int
X2 = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=3;%nov
X3 = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];


figure
scatter3(X1(:,1),X1(:,2),X1(:,3),10,'r.')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),10,'b.')
hold on
scatter3(X3(:,1),X3(:,2),X3(:,3),10,'g.')
hold off

title('Training Distributions')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Novice','Intermediate','Expert')

options = statset('Display','final');
obj1 = fitgmdist(X1,3,'Options',options);
obj2 = fitgmdist(X2,3,'Options',options);
obj3 = fitgmdist(X3,3,'Options',options);

mu1           = mean(X1)
sigma1        = cov(X1)
mu2           = mean(X2)
sigma2        = cov(X2)
mu3           = mean(X3)
sigma3        = cov(X3)


F1 = mvnpdf(X1,mu1,sigma1);
F2 = mvnpdf(X2,mu2,sigma2);
F3 = mvnpdf(X3,mu3,sigma3);

figure
scatter3(X1(:,1),X1(:,2),X1(:,3),8,F1(:));
colormap(cool);
colorbar;
title('Density Estimates: Novice')
xlabel('Position')
ylabel('Speed')
zlabel('Acceleration')
axis([0 30 0 80 0 2*10^5])

figure
scatter3(X2(:,1),X2(:,2),X2(:,3),8,F2(:));
colormap(cool);
colorbar;
title('Density Estimates: Intermediate')
xlabel('Position')
ylabel('Speed')
zlabel('Acceleration')
axis([0 30 0 80 0 2*10^5])

figure
scatter3(X3(:,1),X3(:,2),X3(:,3),8,F3(:));
colormap(cool);
colorbar;
title('Density Estimates: Expert')
xlabel('Position')
ylabel('Speed')
zlabel('Acceleration')
axis([0 30 0 80 0 2*10^5])

%% compute linear  params
jj=1;%nov
dat_nov = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=2;%int
dat_int = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=3;%exp
dat_exp = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Speed(:,1)+All{jj}.Speed(:,2)+All{jj}.Speed(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];

jj=1;%nov
u_nov = abs( All{jj}.Input(:,1)+All{jj}.Input(:,2)+All{jj}.Input(:,3) );
jj=2;%int
u_int = abs( All{jj}.Input(:,1)+All{jj}.Input(:,2)+All{jj}.Input(:,3) );
jj=3;%exp
u_exp = abs( All{jj}.Input(:,1)+All{jj}.Input(:,2)+All{jj}.Input(:,3) );


phi_nov = inv(dat_nov'*dat_nov) * dat_nov' * u_nov;

phi_int = inv(dat_int'*dat_int) * dat_int' * u_int;

phi_exp = inv(dat_exp'*dat_exp) * dat_exp' * u_exp;

%% test linear param classificaiton

myGroups = fliplr([g.gtExp g.flsInt g.flsNov ]);

class = [];

inder = 1;
for gg = 1:length(myGroups) %looping through array of logs
    
    
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)})
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        % ADD SEGMENTS:
        % Extract actual segment data: (left first)
        mySegmentsL = SegScheme{tsk}.dataLogSegs(i, SegType ,L); % Left hand
        
        
        for s = 1:length(mySegmentsL.sgidx) %loop through all segments for that particular grasp   
            
            % Initialize
            tt = [];
            pos = [];
            inputer = [];
            vel = [];
            acc = [];
            
            % selects only the samples of this segment
            iSamples = mySegmentsL.sgidx{s};
            iEnd = iSamples(end);
            iStart = iSamples(1);
            
            tt = DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iEnd, G.Time) ;
            pos = abs( (DataGlb.dataLog{i}( iSamples, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL)) + (DataGlb.dataLog{i}( iSamples, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL)) + (DataGlb.dataLog{i}( iSamples, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL))) ;
            
            %take input as final location of motion
            inputer = abs( (DataGlb.dataLog{i}( iStart, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL))*ones(length(iSamples),1) + (DataGlb.dataLog{i}( iStart, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL))*ones(length(iSamples),1) + (DataGlb.dataLog{i}( iStart, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL))*ones(length(iSamples),1) );
            
            vel = abs( (DataGlb.dataLog{i}( iSamples, G.dxL) - DataGlb.dataLog{i}( iEnd, G.dxL)) + (DataGlb.dataLog{i}( iSamples, G.dyL) - DataGlb.dataLog{i}( iEnd, G.dyL)) + (DataGlb.dataLog{i}( iSamples, G.dzL) - DataGlb.dataLog{i}( iEnd, G.dzL)) );

            TT = mean(diff(tt));
            acc = Calculate_velocity( vel, T, 'holobrodko');
            
            dat = [pos , vel , acc];
            
            err(1) = sum( inputer - dat*phi_nov );
            err(2) = sum( inputer - dat*phi_int );
            err(3) = sum( inputer - dat*phi_exp );
            
            [minn,idx] = min(err);
            
            class(inder,:) = [gg,idx];
            inder = inder + 1;
        end
        
    end
    
    

end
