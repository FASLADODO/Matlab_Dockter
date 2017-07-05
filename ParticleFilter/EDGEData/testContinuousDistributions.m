%test continuous distributions

%http://www.mathworks.com/help/stats/fitgmdist.html

%% first load in this large data structure, if isnt already

HANDLEFT = 1;
HANDRIGHT = 2;

if( ~exist('DataGlb'))
    fprintf('Loading Global Data Structure EdgeDataGlb.mat (huge) ...');
    try
        %load('D:\temp\EdgeDataGlb.mat')
        load('C:\temp\EdgeDataGlb.mat')
    catch
        disp('Could not load local. Loading from M drive...');
        load('M:\Projects\SGP\SGP_DropBoxPort(temp)\Surgery Skills\dataAndAnalysis\Organized Codes\Database\EdgeDataGlb.mat')
    end
    try
        %load('D:\temp\EDGE_Segments.mat')
        load('C:\temp\EDGE_Segments.mat')
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

%delete
clear All Segment


%find minimum number of segments 
totalsegs = [];
for gg = 1:length(myGroups) %looping through array of logs
    segind = 0;
    myGroups(gg)
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)})
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        % ADD SEGMENTS:
        % Extract actual segment data: (left first)
        mySegmentsL = SegScheme{tsk}.dataLogSegs(i, SegType ,L); % Left hand
        
        for s = 1:length(mySegmentsL.sgidx) %loop through all segments for that particular grasp 
            segind = segind + 1;
            
        end
    end
    
    totalsegs(gg) = segind;
end

totalsegs
limitseg = min(totalsegs)


trajlen = [];

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
    
    meanlen = [];
    
    trj_index = 1;
    segind = 0;
    
    % Sum up all logs specific group (int, nov, exp)
    for ii = 1:length(DataGlb.grp.all{tsk}.Idx{myGroups(gg)})
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        % ADD SEGMENTS:
        % Extract actual segment data: (left first)
        mySegmentsL = SegScheme{tsk}.dataLogSegs(i, SegType ,L); % Left hand
        
        for s = 1:length(mySegmentsL.sgidx) %loop through all segments for that particular grasp   
            
            if(segind < limitseg)
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
            %incrment counter
            segind = segind + 1;
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

    trajlen{gg} = meanlen;
end



mean(trajlen{1})
mean(trajlen{2})
mean(trajlen{3})

%% 2D probabalities 

xyz = 3; %z axis
az = 20;
el = 30;

limz = [-30 30 -100 100];

jj = 1;
X1 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];
jj = 2;
X2 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];
jj = 3;
X3 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];

%get gaussian mixture model
GMModel1 = fitgmdist(X1,1);
GMModel2 = fitgmdist(X2,1);
GMModel3 = fitgmdist(X3,1);


GMModel1.Sigma
GMModel2.Sigma
GMModel3.Sigma


NP = 50;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);

Xon = [XG1,XG2];


F1 = pdf(GMModel1,Xon);
F2 = pdf(GMModel2,Xon);
F3 = pdf(GMModel3,Xon);


F1 = F1 ./ max(F1);
F2 = F2 ./ max(F2);
F3 = F3 ./ max(F3);

figure
scatter(X1(:,1),X1(:,2),2);
hold on
m1 = Surface3Dalt(Xon(:,1),Xon(:,2),F1,limz,NP);
set(m1,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF Z (exp)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);

figure
scatter(X2(:,1),X2(:,2),2);
hold on
m2 = Surface3Dalt(Xon(:,1),Xon(:,2),F2,limz,NP);
set(m2,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF Z (int)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);


figure
scatter(X3(:,1),X3(:,2),2);
hold on
m3 = Surface3Dalt(Xon(:,1),Xon(:,2),F3,limz,NP);
set(m3,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF Z (nov)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);


%Now see max diffs
NP = 1000;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);

Xon = [XG1,XG2];

F1 = pdf(GMModel1,Xon);
F2 = pdf(GMModel2,Xon);
F3 = pdf(GMModel3,Xon);

F1 = F1 ./ max(F1);
F2 = F2 ./ max(F2);
F3 = F3 ./ max(F3);

ERR1 = abs(F1 - F2);
ERR2 = abs(F1 - F3);
ERR3 = abs(F3 - F2);

seperatz(xyz,:) = [max(ERR1), max(ERR2), max(ERR3)]

AllDist{xyz}.Class{1} = GMModel1;
AllDist{xyz}.Class{2} = GMModel2;
AllDist{xyz}.Class{3} = GMModel3;


%% 2D probabalities 

xyz = 2; %y axis
az = 20;
el = 30;

limz = [-30 30 -100 100];

jj = 1;
X1 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];
jj = 2;
X2 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];
jj = 3;
X3 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];

%get gaussian mixture model
GMModel1 = fitgmdist(X1,1);
GMModel2 = fitgmdist(X2,1);
GMModel3 = fitgmdist(X3,1);


NP = 50;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);

Xon = [XG1,XG2];

F1 = pdf(GMModel1,Xon);
F2 = pdf(GMModel2,Xon);
F3 = pdf(GMModel3,Xon);

F1 = F1 ./ max(F1);
F2 = F2 ./ max(F2);
F3 = F3 ./ max(F3);

figure
scatter(X1(:,1),X1(:,2),2);
hold on
m1 = Surface3Dalt(Xon(:,1),Xon(:,2),F1,limz,NP);
set(m1,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF Y (exp)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);

figure
scatter(X2(:,1),X2(:,2),2);
hold on
m2 = Surface3Dalt(Xon(:,1),Xon(:,2),F2,limz,NP);
set(m2,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF Y (int)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);


figure
scatter(X3(:,1),X3(:,2),2);
hold on
m3 = Surface3Dalt(Xon(:,1),Xon(:,2),F3,limz,NP);
set(m3,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF Y (nov)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);

%Now see max diffs
NP = 1000;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);

Xon = [XG1,XG2];

F1 = pdf(GMModel1,Xon);
F2 = pdf(GMModel2,Xon);
F3 = pdf(GMModel3,Xon);

F1 = F1 ./ max(F1);
F2 = F2 ./ max(F2);
F3 = F3 ./ max(F3);

ERR1 = abs(F1 - F2);
ERR2 = abs(F1 - F3);
ERR3 = abs(F3 - F2);

seperatz(xyz,:) = [max(ERR1), max(ERR2), max(ERR3)];

AllDist{xyz}.Class{1} = GMModel1;
AllDist{xyz}.Class{2} = GMModel2;
AllDist{xyz}.Class{3} = GMModel3;



%% 2D probabalities 

xyz = 1; %x axis
az = 20;
el = 30;

limz = [-30 30 -100 100];

jj = 1;
X1 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];
jj = 2;
X2 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];
jj = 3;
X3 = [All{jj}.Vel(:,xyz), All{jj}.Acc(:,xyz)];

%get gaussian mixture model
GMModel1 = fitgmdist(X1,1);
GMModel2 = fitgmdist(X2,1);
GMModel3 = fitgmdist(X3,1);


NP = 50;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);

Xon = [XG1,XG2];

F1 = pdf(GMModel1,Xon);
F2 = pdf(GMModel2,Xon);
F3 = pdf(GMModel3,Xon);

F1 = F1 ./ max(F1);
F2 = F2 ./ max(F2);
F3 = F3 ./ max(F3);

figure
scatter(X1(:,1),X1(:,2),2);
hold on
m1 = Surface3Dalt(Xon(:,1),Xon(:,2),F1,limz,NP);
set(m1,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF X (exp)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);

figure
scatter(X2(:,1),X2(:,2),2);
hold on
m2 = Surface3Dalt(Xon(:,1),Xon(:,2),F2,limz,NP);
set(m2,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF X (int)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);


figure
scatter(X3(:,1),X3(:,2),2);
hold on
m3 = Surface3Dalt(Xon(:,1),Xon(:,2),F3,limz,NP);
set(m3,'facecolor','none')
hold off
colormap(cool);
colorbar;
title('Scatter Plot and Fitted PDF X (nov)')
xlabel('vel')
ylabel('acc')
zlabel('probs')
axis([limz 0 1 ])
view(az, el);

%Now see max diffs
NP = 1000;
xo1 = linspace(limz(1),limz(2),NP);
xo2 = linspace(limz(3),limz(4),NP);
[XGt1,XGt2] = meshgrid(xo1,xo2);
XG1 = reshape(XGt1,[],1);
XG2 = reshape(XGt2,[],1);

Xon = [XG1,XG2];

F1 = pdf(GMModel1,Xon);
F2 = pdf(GMModel2,Xon);
F3 = pdf(GMModel3,Xon);

F1 = F1 ./ max(F1);
F2 = F2 ./ max(F2);
F3 = F3 ./ max(F3);

ERR1 = abs(F1 - F2);
ERR2 = abs(F1 - F3);
ERR3 = abs(F3 - F2);

seperatz(xyz,:) = [max(ERR1), max(ERR2), max(ERR3)];

AllDist{xyz}.Class{1} = GMModel1;
AllDist{xyz}.Class{2} = GMModel2;
AllDist{xyz}.Class{3} = GMModel3;


%% Try classifying with probs

xyz = 3; %z

thresh = 0.05;

for gg = 1:length(Segment)
   for tt = 1:length(Segment{gg}.Trajectory)
       for ii = 1:length(Segment{gg}.Trajectory{tt}.Vel)
           xon = [Segment{gg}.Trajectory{tt}.Vel(ii,xyz), Segment{gg}.Trajectory{tt}.Acc(ii,xyz)];
           
           probz = [];
           for cc = 1:length(AllDist{xyz}.Class)
               probz(cc) = pdf(AllDist{xyz}.Class{cc},xon);
           end
           err1_temp = abs(probz(1) - probz(2));
           err2_temp = abs(probz(1) - probz(3));
           err3_temp = abs(probz(2) - probz(3));
           
       end
   end
end


%% 3D probabalities

X1 = All{1}.Acc;
X2 = All{2}.Acc;
X3 = All{3}.Acc;

%get gaussian mixture model
GMModel1 = fitgmdist(X1,1);
GMModel2 = fitgmdist(X2,1);
GMModel3 = fitgmdist(X3,1);

F1 = pdf(GMModel1,X1);
F2 = pdf(GMModel2,X2);
F3 = pdf(GMModel3,X3);


F1 = F1 ./ max(F1);
F2 = F2 ./ max(F2);
F3 = F3 ./ max(F3);

figure
hs1 = scatter3(X1(:,1),X1(:,2),X1(:,3),2,F1(:));
colormap(cool);
colorbar;
alpha(0.3);
title('Density Estimates Expert')
xlabel('ddx')
ylabel('ddy')
zlabel('ddz')


figure
hs2 = scatter3(X2(:,1),X2(:,2),X2(:,3),2,F2(:));
colormap(cool);
colorbar;
alpha(0.3);
title('Density Estimates Intermdiate')
xlabel('ddx')
ylabel('ddy')
zlabel('ddz')

figure
hs3 = scatter3(X3(:,1),X3(:,2),X3(:,3),2,F3(:));
colormap(cool);
colorbar;
alpha(0.3);
title('Density Estimates Novice')
xlabel('ddx')
ylabel('ddy')
zlabel('ddz')
