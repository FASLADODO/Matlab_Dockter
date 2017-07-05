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

HANDLEFT = 1;
HANDRIGHT = 2;

if( ~exist('DataGlb'))
    fprintf('Loading Global Data Structure EdgeDataGlb.mat (huge) ...');
    try
        load('D:\temp\EdgeDataGlb.mat')
        %load('C:\temp\EdgeDataGlb.mat')
    catch
        disp('Could not load local. Loading from M drive...');
        load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EdgeDataGlb.mat')
    end
    try
        load('D:\temp\EDGE_Segments.mat')
        %load('C:\temp\EDGE_Segments.mat')
    catch
        disp('Could not load local. Loading from M drive...');
        load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EDGE_Segments.mat')
    end
    fprintf(' DONE.\n');
end

tsk = 1; % Peg Transfer task
myGroups = [g.gtExp g.flsInt g.flsNov ]; % select only three skill groups
classGroup = {'expert';'intermediate';'novice'};


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

meanlen = [];

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
    
    trj_index = 1;
    
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
            
            
            
            %%%%%%%%%%%%%%%%%%%get time and pos and velocity
            tempt = DataGlb.dataLog{i}( iSamples, G.Time) - DataGlb.dataLog{i}( iStart, G.Time);
            tempx = DataGlb.dataLog{i}( iSamples, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL);
            tempy = DataGlb.dataLog{i}( iSamples, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL);
            tempz = DataGlb.dataLog{i}( iSamples, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL);
            
            %take input as distance to final position
            tempinx = (DataGlb.dataLog{i}( iStart, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL))*ones(length(iSamples),1); 
            tempiny = (DataGlb.dataLog{i}( iStart, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL))*ones(length(iSamples),1); 
            tempinz = (DataGlb.dataLog{i}( iStart, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL))*ones(length(iSamples),1);
            
            tempinx_pulse = zeros(length(iSamples),1); 
            tempiny_pulse = zeros(length(iSamples),1); 
            tempinz_pulse = zeros(length(iSamples),1);
            tempinx_pulse(1,:) = DataGlb.dataLog{i}( iStart, G.xL) - DataGlb.dataLog{i}( iEnd, G.xL); 
            tempiny_pulse(1,:) = DataGlb.dataLog{i}( iStart, G.yL) - DataGlb.dataLog{i}( iEnd, G.yL); 
            tempinz_pulse(1,:) = DataGlb.dataLog{i}( iStart, G.zL) - DataGlb.dataLog{i}( iEnd, G.zL);
            
            tempdx = DataGlb.dataLog{i}( iSamples, G.dxL);
            tempdy = DataGlb.dataLog{i}( iSamples, G.dyL);
            tempdz = DataGlb.dataLog{i}( iSamples, G.dzL);
            
            T = mean(diff(tempt));
            tempddx = Calculate_velocity( tempdx, T, 'holobrodko') ;
            tempddy = Calculate_velocity( tempdy, T, 'holobrodko') ;
            tempddz = Calculate_velocity( tempdz, T, 'holobrodko') ;
            
            tempgrp = repmat( DataGlb.grp.tag(myGroups(gg)), [length(iSamples),1]);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%Store individual grasps
            meanlen = [meanlen, length(tempx)];
            Segment{gg}.Trajectory{trj_index}.Input_All = [tempinx,tempiny,tempinz];
            Segment{gg}.Trajectory{trj_index}.Input = sqrt(tempinx_pulse.^2 + tempiny_pulse.^2 + tempinz_pulse.^2) ;
            Segment{gg}.Trajectory{trj_index}.State_All = [tempx,tempy,tempz,tempdx,tempdy,tempdz,tempddx,tempddy,tempddz];
            Segment{gg}.Trajectory{trj_index}.State = [abs(tempx.^2 + tempy.^2 + tempz.^2),abs(tempdx.^2 + tempdy.^2 + tempdz.^2),abs(tempddx.^2 + tempddy.^2 + tempddz.^2)];
            Segment{gg}.Trajectory{trj_index}.Pos = [tempx,tempy,tempz];
            Segment{gg}.Trajectory{trj_index}.Vel = [tempdx,tempdy,tempdz];
            Segment{gg}.Trajectory{trj_index}.Acc = [tempddx,tempddy,tempddz];
            Segment{gg}.Trajectory{trj_index}.Grp = tempgrp;
            Segment{gg}.Trajectory{trj_index}.Time = tempt;
            trj_index = trj_index + 1;
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%store all into struct
            t = [ t ; tempt] ;
            xleft = [ xleft ; tempx] ;
            yleft = [ yleft ; tempy] ;
            zleft = [ zleft ; tempz] ;

            input_x = [input_x; tempinx];
            input_y = [input_y; tempiny];
            input_z = [input_z; tempinz];

            dxleft = [ dxleft ; tempdx ] ;
            dyleft = [ dyleft ; tempdy ] ;
            dzleft = [ dzleft ; tempdz ] ;

            ddxleft = [ ddxleft ; tempddx ] ;
            ddyleft = [ ddyleft ; tempddy ] ;
            ddzleft = [ ddzleft ; tempddz ] ;



            grpL = [ grpL; tempgrp ];


        end
        
    end
    
    InputL = [input_x, input_y, input_z ];
    PosL = [xleft, yleft, zleft ];
    VelL = [dxleft, dyleft, dzleft ];
    AccL = [ddxleft, ddyleft, ddzleft ];
    
    All{gg}.Input = InputL;
    All{gg}.Vel = VelL;
    All{gg}.Acc = AccL;
    All{gg}.Pos  =  PosL;
    All{gg}.Grp   = grpL ;


end

mean(meanlen)
min(meanlen)
max(meanlen)
std(meanlen)

%% Play with this stuff:

%http://www.mathworks.com/help/stats/continuous-distributions.html

%% Play wif dat data

[nn,order_s] = size(Segment{1}.Trajectory{1}.State_All)
[nn,order_in] = size(Segment{1}.Trajectory{1}.Input_All)

%get parameters (should be 9x3 parameters) 3 inputs, 9 states
for gg = 1:length(Segment)
    Param_General{gg} = [];
    for tt = 1:length(Segment{gg}.Trajectory)
        %input_temp = sqrt(sum(Segment{gg}.Trajectory{tt}.Input.^2,2)); 
        paramtemp = pinv(Segment{gg}.Trajectory{tt}.State_All)*Segment{gg}.Trajectory{tt}.Input_All;
        Param_General{gg} = [Param_General{gg}, reshape(paramtemp,3*order_s,1) ];
    end
    temp_p_Shape = mean(Param_General{gg},2);
    Param_AVG{gg} = reshape(temp_p_Shape,order_s,order_in);
end

%try classifying with big ol parameters
Error_Simple = [];
for gg = 1:length(Segment)
    indyr = 1;
    for tt = 1:length(Segment{gg}.Trajectory)
        for cc = 1:length(Segment)
            %input_temp = sqrt(sum(Segment{gg}.Trajectory{tt}.Input.^2,2));
            Error_Simple{gg}(indyr,cc) = sum(sum(abs(Segment{gg}.Trajectory{tt}.Input_All - Segment{gg}.Trajectory{tt}.State_All*Param_AVG{cc})));
        
        end
        [val,idx] = min(Error_Simple{gg}(indyr,:));
        Classificate{gg}(indyr,:) = [gg,idx];
        indyr = indyr + 1;
    end
    
    %get accuracy
    acctemp = Classificate{gg}(:,1) == Classificate{gg}(:,2);
    Accuracy{gg} = sum(acctemp)/length(acctemp);
end


meas = [];
skillLevel = [];
indall = 1;

%now trying a ratio of time to overall path length?
for gg = 1:length(Segment)
    indyr = 1;
    for tt = 1:length(Segment{gg}.Trajectory)
        ratio = abs( Segment{gg}.Trajectory{tt}.Time(end) / Segment{gg}.Trajectory{tt}.Input(1) );
        %ratio2 = abs( max(mean(Segment{gg}.Trajectory{tt}.Vel )) / Segment{gg}.Trajectory{tt}.Input(1) );
        ratio2 = mean( abs( Segment{gg}.Trajectory{tt}.Acc ) );
        %variance in velocity/acc
        ratio3 = var( Segment{gg}.Trajectory{tt}.Vel );
        ratio4 = var( Segment{gg}.Trajectory{tt}.Acc );
        %sign changes in velocity
        ratio5 = signChanges( Segment{gg}.Trajectory{tt}.Vel );
        Ratio_All{gg}(indyr) = ratio;
        Ratio_All_2{gg}(indyr,:) = ratio2;
        Ratio_All_5{gg}(indyr,:) = ratio5;
        indyr = indyr + 1;
        
        % For LDA
        meas(indall,:) = [ratio,ratio5];
        skillLevel{indall,1} = classGroup{gg};
        indall = indall + 1;
    end
    
%     Ratio_Mean(gg) = mean(Ratio_All{gg});
%     Ratio_Std(gg) = std(Ratio_All{gg});
%     Ratio_Mean_2(gg) = mean(Ratio_All_2{gg});
%     Ratio_Std_2(gg) = std(Ratio_All_2{gg});
     Ratio_Mean_5(gg,:) = mean(Ratio_All_5{gg});
end

Ratio_Mean_5

%get classifier
quadclass = fitcdiscr(meas,skillLevel,'discrimType','quadratic');

%test classify
for ii = 1:length(meas);
    meas_online = meas(ii,:);
    class_online = predict(quadclass,meas_online);
    classinyo{ii,1} = skillLevel{ii,1};
    classinyo{ii,2} = class_online;
    
    acc(ii) = strcmp(classinyo{ii,1},classinyo{ii,2});
end

sum(acc)/length(acc)

%% Plotting QDA


PL = meas(:,2);
PW = meas(:,1);

h1 = gscatter(PL,PW,skillLevel,'krb','ov^',[],'off');
h1(1).LineWidth = 2;
h1(2).LineWidth = 2;
h1(3).LineWidth = 2;
legend('Nov','Int','Exp','Location','best')
hold on

X = [PL,PW];
cqs = fitcdiscr(X,skillLevel,...
    'DiscrimType','quadratic');

% Now, retrieve the coefficients for the quadratic boundary between the
% second and third classes (versicolor and virginica).
K = cqs.Coeffs(2,3).Const;
L = cqs.Coeffs(2,3).Linear;
Q = cqs.Coeffs(2,3).Quadratic;

% Plot the curve K + [x1,x2]*L + [x1,x2]*Q*[x1,x2]' = 0.
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;
h2 = ezplot(f,[0 20 0 7]);
h2.Color = 'r';
h2.LineWidth = 2;

% Now, retrieve the coefficients for the quadratic boundary between the
% first and second classes (setosa and versicolor).
K = cqs.Coeffs(1,2).Const;
L = cqs.Coeffs(1,2).Linear;
Q = cqs.Coeffs(1,2).Quadratic;

% Plot the curve K + [x1,x2]*L + [x1,x2]*Q*[x1,x2]'=0:
f = @(x1,x2) K + L(1)*x1 + L(2)*x2 + Q(1,1)*x1.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + Q(2,2)*x2.^2;
h3 = ezplot(f,[0 20 0 7]); % Plot the relevant portion of the curve.
h3.Color = 'k';
h3.LineWidth = 2;
axis([0 20 0 7])
xlabel('path')
ylabel('vel var')
title('{\bf Quadratic Classification with EDGE Training Data}')
hold off


%% Plot ratio 1 vs ratio 3

figure
gg = 1;
scatter(Ratio_All_5{gg}(:,3),Ratio_All{gg},10,'r.')
hold on
gg = 2;
scatter(Ratio_All_5{gg}(:,3),Ratio_All{gg},10,'b.')
hold on
gg = 3;
scatter(Ratio_All_5{gg}(:,3),Ratio_All{gg},10,'g.')
hold off

title('Ratio Distributions')
xlabel('Class')
ylabel('Ratio')
legend('Expert','Intermediate','Novice')


%% Plot ratio 1 to find thresholds

figure
gg = 1;
scatter(ones(length(Ratio_All{gg}),1)*gg,Ratio_All{gg},10,'r.')
hold on
gg = 2;
scatter(ones(length(Ratio_All{gg}),1)*gg,Ratio_All{gg},10,'b.')
hold on
gg = 3;
scatter(ones(length(Ratio_All{gg}),1)*gg,Ratio_All{gg},10,'g.')
hold off

title('Ratio Distributions')
xlabel('Class')
ylabel('Ratio')
legend('Expert','Intermediate','Novice')

%% Plot ratio 5 to find thresholds

figure
gg = 1;
scatter(ones(length(Ratio_All_5{gg}),1)*gg,Ratio_All_5{gg}(:,3),10,'r.')
hold on
gg = 2;
scatter(ones(length(Ratio_All_5{gg}),1)*gg,Ratio_All_5{gg}(:,3),10,'b.')
hold on
gg = 3;
scatter(ones(length(Ratio_All_5{gg}),1)*gg,Ratio_All_5{gg}(:,3),10,'g.')
hold off

title('Ratio Distributions')
xlabel('Class')
ylabel('Ratio')
legend('Expert','Intermediate','Novice')


%% Plot pos with bounding spheres

spherecolor = [1 1 0];

alph = 0.99;

%expert
figure
jj = 1; %exp
dat1 = All{jj}.Pos;
[radius, center] = bounding_sphere(dat1,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat1(:,1),dat1(:,2),dat1(:,3),'r.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Position Sphere Expert')
xlabel('x')
ylabel('y')
zlabel('z')
legend('pos','bounding')

%store that
SphereB{jj}.Pos = radius;
CenterB{jj}.Pos = center;

%check if actually 95%
poo = [];
daton = dat1; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy

%Intermediate
figure
jj = 2; %int
dat2 = All{jj}.Pos;
[radius, center] = bounding_sphere(dat2,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat2(:,1),dat2(:,2),dat2(:,3),'b.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Position Sphere Intermediate')
xlabel('x')
ylabel('y')
zlabel('z')
legend('pos','bounding')

%store that
SphereB{jj}.Pos = radius;
CenterB{jj}.Pos = center;

%check if actually 95%
poo = [];
daton = dat2; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy

%novice
figure
jj = 3; %nov
dat3 = All{jj}.Pos;
[radius, center] = bounding_sphere(dat3,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat3(:,1),dat3(:,2),dat3(:,3),'g.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Position Sphere Novice')
xlabel('x')
ylabel('y')
zlabel('z')
legend('pos','bounding')

%store that
SphereB{jj}.Pos = radius;
CenterB{jj}.Pos = center;

%check if actually 95%
poo = [];
daton = dat3; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy

%% Plot Velocity with bounding spheres

spherecolor = [1 1 0];
alph = 0.99;

%expert
figure
jj = 1; %exp
dat1 = All{jj}.Vel;
[radius, center] = bounding_sphere(dat1,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat1(:,1),dat1(:,2),dat1(:,3),'r.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Velocity Sphere Expert')
xlabel('dx')
ylabel('dy')
zlabel('dz')
legend('vel','bounding')

%store that
SphereB{jj}.Vel = radius;
CenterB{jj}.Vel = center;

%check if actually 95%
poo = [];
daton = dat1; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy

%Intermediate
figure
jj = 2; %int
dat2 = All{jj}.Vel;
[radius, center] = bounding_sphere(dat2,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat2(:,1),dat2(:,2),dat2(:,3),'b.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Velocity Sphere Intermediate')
xlabel('dx')
ylabel('dy')
zlabel('dz')
legend('vel','bounding')

%store that
SphereB{jj}.Vel = radius;
CenterB{jj}.Vel = center;

%check if actually 95%
poo = [];
daton = dat2; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy

%novice
figure
jj = 3; %nov
dat3 = All{jj}.Vel;
[radius, center] = bounding_sphere(dat3,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat3(:,1),dat3(:,2),dat3(:,3),'g.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Velocity Sphere Novice')
xlabel('dx')
ylabel('dy')
zlabel('dz')
legend('vel','bounding')

%store that
SphereB{jj}.Vel = radius;
CenterB{jj}.Vel = center;

%check if actually 95%
poo = [];
daton = dat3; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy


%% Plot Acceleration with bounding spheres

spherecolor = [1 1 0];
alph = 0.99;

%expert
figure
jj = 1; %exp
dat1 = All{jj}.Acc;
[radius, center] = bounding_sphere(dat1,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat1(:,1),dat1(:,2),dat1(:,3),'r.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Velocity Sphere Expert')
xlabel('dx')
ylabel('dy')
zlabel('dz')
legend('Acceleration','bounding')

%store that
SphereB{jj}.Acc = radius;
CenterB{jj}.Acc = center;
RadiusAllAcc(jj) = radius;

%check if actually 95%
poo = [];
daton = dat1; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy

%Intermediate
figure
jj = 2; %int
dat2 = All{jj}.Acc;
[radius, center] = bounding_sphere(dat2,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat2(:,1),dat2(:,2),dat2(:,3),'b.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Acceleration Sphere Intermediate')
xlabel('dx')
ylabel('dy')
zlabel('dz')
legend('vel','bounding')

%store that
SphereB{jj}.Acc = radius;
CenterB{jj}.Acc = center;
RadiusAllAcc(jj) = radius;

%check if actually 95%
poo = [];
daton = dat2; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy

%novice
figure
jj = 3; %nov
dat3 = All{jj}.Acc;
[radius, center] = bounding_sphere(dat3,alph)
[X,Y,Z] = sphereObject(radius,center);
scatter3(dat3(:,1),dat3(:,2),dat3(:,3),'g.')
hold on
m=mesh(X,Y,Z);
set(m,'facecolor','none')
colormap(spherecolor)
hold off
title('Acceleration Sphere Novice')
xlabel('dx')
ylabel('dy')
zlabel('dz')
legend('vel','bounding')

%store that
SphereB{jj}.Acc = radius;
CenterB{jj}.Acc = center;
RadiusAllAcc(jj) = radius;

%check if actually 95%
poo = [];
daton = dat3; %get data
datc = daton - center(ones(size(daton,1),1),:); %centerthatshit
mag = sqrt(sum(datc.^2, 2)); %distance
radz = ones(length(datc),1)*radius; %vector of radius for comparing
poo = mag <= radz; %compare it
acc = sum(poo)/length(poo) %accuracy


%% Try classifying with sphere bounds

alph = 0.99;

indd = 1;

for gg = 1:length(Segment)
    
    for tt = 1:length(Segment{gg}.Trajectory)
        datp = Segment{gg}.Trajectory{tt}.Pos;
        datv = Segment{gg}.Trajectory{tt}.Vel;
        data = Segment{gg}.Trajectory{tt}.Acc;
        [radp, cenp] = bounding_sphere(datp,alph);
        [radv, cenv] = bounding_sphere(datv,alph);
        [rada, cena] = bounding_sphere(data,alph);
        AvgSphere{gg}.Pos(indd) = radp;
        AvgSphere{gg}.Vel(indd) = radv;
        AvgSphere{gg}.Acc(indd) = rada;
        indd = indd + 1;
        
        err = RadiusAllAcc - rada;
        [val,coridx] = min(abs(err));
        classin(indd,:) = [gg,coridx];
    end
end

error = classin(:,1) == classin(:,2);
accuracy = sum(error)/length(error)



%% Try x y z phase portraits seperate

xyz = 3;

figure
jj=1;%exp
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'r.');
hold on
jj=2;%int
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'b.');
hold on
jj=3;%nov
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'g.');
hold off

title('Z axis phase portrait')
xlabel('Vel')
ylabel('Acc')
legend('Expert','Intermediate','Novice')


xyz = 2;

figure
jj=1;%exp
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'r.');
hold on
jj=2;%int
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'b.');
hold on
jj=3;%nov
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'g.');
hold off

title('Y axis phase portrait')
xlabel('Vel')
ylabel('Acc')
legend('Expert','Intermediate','Novice')

xyz = 1;

figure
jj=1;%exp
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'r.');
hold on
jj=2;%int
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'b.');
hold on
jj=3;%nov
scatter(All{jj}.Vel(:,xyz),All{jj}.Acc(:,xyz),10,'g.');
hold off

title('X axis phase portrait')
xlabel('Vel')
ylabel('Acc')
legend('Expert','Intermediate','Novice')



%% Plot pos and vel and acc on different plots


figure
jj=1;%exp
scatter3(All{jj}.Pos(:,1),All{jj}.Pos(:,2),All{jj}.Pos(:,3),10,'r.');
hold on
jj=2;%int
scatter3(All{jj}.Pos(:,1),All{jj}.Pos(:,2),All{jj}.Pos(:,3),10,'b.');
hold on
jj=3;%nov
scatter3(All{jj}.Pos(:,1),All{jj}.Pos(:,2),All{jj}.Pos(:,3),10,'g.');
hold off

title('Position 3D')
xlabel('X')
ylabel('Y')
zlabel('Z')
legend('Expert','Intermediate','Novice')

figure
jj=1;%exp
scatter3(All{jj}.Vel(:,1),All{jj}.Vel(:,2),All{jj}.Vel(:,3),10,'r.');
hold on
jj=2;%int
scatter3(All{jj}.Vel(:,1),All{jj}.Vel(:,2),All{jj}.Vel(:,3),10,'b.');
hold on
jj=3;%now
scatter3(All{jj}.Vel(:,1),All{jj}.Vel(:,2),All{jj}.Vel(:,3),10,'g.');
hold off

title('Velocity 3D')
xlabel('dX')
ylabel('dY')
zlabel('dZ')
legend('Expert','Intermediate','Novice')

figure
jj=1;%exp
scatter3(All{jj}.Acc(:,1),All{jj}.Acc(:,2),All{jj}.Acc(:,3),10,'r.');
hold on
jj=2;%int
scatter3(All{jj}.Acc(:,1),All{jj}.Acc(:,2),All{jj}.Acc(:,3),10,'b.');
hold on
jj=3;%nov
scatter3(All{jj}.Acc(:,1),All{jj}.Acc(:,2),All{jj}.Acc(:,3),10,'g.');
hold off

title('Acceleration 3D')
xlabel('ddX')
ylabel('ddY')
zlabel('ddZ')
legend('Expert','Intermediate','Novice')



%% Plot pos,vel,acc on same plot for all classes

jj=1;%exp
dat_nov = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Vel(:,1)+All{jj}.Vel(:,2)+All{jj}.Vel(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=2;%int
dat_int = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Vel(:,1)+All{jj}.Vel(:,2)+All{jj}.Vel(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];
jj=3;%nov
dat_exp = [abs( All{jj}.Pos(:,1)+All{jj}.Pos(:,2)+All{jj}.Pos(:,3) ) , abs( All{jj}.Vel(:,1)+All{jj}.Vel(:,2)+All{jj}.Vel(:,3) ), abs( All{jj}.Acc(:,1)+All{jj}.Acc(:,2)+All{jj}.Acc(:,3) ) ];

dat = [dat_nov; dat_int; dat_exp];

figure
scatter3(dat_nov(:,1),dat_nov(:,2),dat_nov(:,3),10,'r.')
hold on
scatter3(dat_int(:,1),dat_int(:,2),dat_int(:,3),10,'b.')
hold on
scatter3(dat_exp(:,1),dat_exp(:,2),dat_exp(:,3),10,'g.')
hold off

title('Training Distributions')
xlabel('Position')
ylabel('Velcocity')
zlabel('Acceleration')
legend('Expert','Intermediate','Novice')




%% plot each class on seperate plots

sampz = 0.1;


jj=1;%exp
len = length(Segment{jj}.Trajectory);
subsam = randperm(len);
SS = round(len*sampz);
dat_exp_s = [];
for ii = 1:SS
    strct = Segment{jj}.Trajectory{subsam(ii)};
    dat_exp_s = [dat_exp_s ; abs( strct.Pos(:,1)+strct.Pos(:,2)+strct.Pos(:,3) ) , abs( strct.Vel(:,1)+strct.Vel(:,2)+strct.Vel(:,3) ), abs( strct.Acc(:,1)+strct.Acc(:,2)+strct.Acc(:,3) ) ];
end

jj=2;%int
len = length(Segment{jj}.Trajectory);
subsam = randperm(len);
SS = round(len*sampz);
dat_int_s = [];
for ii = 1:SS
    strct = Segment{jj}.Trajectory{subsam(ii)};
    dat_int_s = [dat_int_s ; abs( strct.Pos(:,1)+strct.Pos(:,2)+strct.Pos(:,3) ) , abs( strct.Vel(:,1)+strct.Vel(:,2)+strct.Vel(:,3) ), abs( strct.Acc(:,1)+strct.Acc(:,2)+strct.Acc(:,3) ) ];
end

jj=3;%nov
len = length(Segment{jj}.Trajectory);
subsam = randperm(len);
SS = round(len*sampz);
dat_nov_s = [];
for ii = 1:SS
    strct = Segment{jj}.Trajectory{subsam(ii)};
    dat_nov_s = [dat_nov_s ; abs( strct.Pos(:,1)+strct.Pos(:,2)+strct.Pos(:,3) ) , abs( strct.Vel(:,1)+strct.Vel(:,2)+strct.Vel(:,3) ), abs( strct.Acc(:,1)+strct.Acc(:,2)+strct.Acc(:,3) ) ];
end

figure
scatter3(dat_nov_s(:,1),dat_nov_s(:,2),dat_nov_s(:,3),10,'r.')
title('Training Novice')
xlabel('Position')
ylabel('Velcocity')
zlabel('Acceleration')
axis([0,25,0,60,0,200])

figure
scatter3(dat_int_s(:,1),dat_int_s(:,2),dat_int_s(:,3),10,'b.')
title('Training Intermediate')
xlabel('Position')
ylabel('Velcocity')
zlabel('Acceleration')
axis([0,25,0,60,0,200])

figure
scatter3(dat_exp_s(:,1),dat_exp_s(:,2),dat_exp_s(:,3),10,'g.')
title('Training Expert')
xlabel('Position')
ylabel('Velcocity')
zlabel('Acceleration')
axis([0,25,0,60,0,200])





