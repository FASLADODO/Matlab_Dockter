%Test RF with EDGE data, using only the sub dimension from reliefrbf

% first load in this large data structure, if isnt already

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
fsize = 14;

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
    for ii = 1:allsurgeons(gg)
        
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

%% Create dummy Data and labels

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

key.strings = {'dx','dy','dz','ddx','ddy','ddz','$\ddot{x}$','$\ddot{y}$','$\ddot{z}$','velmag','velalpha','velbeta','accmag','accalpha','accbeta','$\theta$','$\dot{\theta}$'};
%We use ddot instead of dddot beacue it doesnt exists
% key.strings = {'$\ddot{x}$','$\ddot{y}$','$\ddot{z}$','$\theta$','$\dot{\theta}$'};

%which states we want to keep
StrLabels = {'Novice';'Intermediate';'Expert'};
grplist = [1,3];
mapping = [0,-1,1];

%clear this monster
DataAll = [];
LabelsAll = [];

%loop through all groups trials, surgeons, hands, segments
for gg = grplist %skill levels
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};

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

                DataAll = [DataAll; dtemp];
                LabelsAll = [LabelsAll; ones(ng,1)*mapping(gg)];
                
                SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.AllStates = dtemp;
            end
        end
    end
end

%test columns
tc1 = key.col.grip ;
tc2 = key.col.dddx;
tc3 = key.col.dddz;

figure
gscatter3(DataAll(:,tc1),DataAll(:,tc2),DataAll(:,tc3),LabelsAll,{'c','r'},{'o','+'})
% title('phase portrait train')
xlabel(key.strings(tc1),'FontSize',fsize,'Interpreter','latex');
ylabel(key.strings(tc2),'FontSize',fsize,'Interpreter','latex');
zlabel(key.strings(tc3),'FontSize',fsize,'Interpreter','latex');
leg = findobj(gcf,'Tag','legend');
set(leg,'FontSize',fsize);


%% leave one user out 

% get data
runs = allsurgeons(1); %the number of novice runs

%for prettiness
colormapnew = flipud(cool);
trees = 20;

%state and input functions (xdot + x + x^2
coluse = [key.col.grip, key.col.dddx, key.col.dddz];
% coluse = [key.col.accmag, key.col.velmag, key.col.grip]; %pearsons
mapping = [0,-1,1];

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
ScoreTime = [];
DataStore = [];


for oo = 1:runs %length(SegData.Donor)
    
    fprintf('Epoch %d of %d \n',oo,runs)
    
    %clear our matrices
    Xtrain = [];
    Labelstrain = [];
    Xtest = [];
    Labelstest = [];
    Ptrain = []; %parameters
    
    %stash just the training data
    for gg = grplist %skill levels
        numtrial = length(SegData{gg}.Trial);
        leaveout = mod(oo-1,numtrial)+1; %mod 1 indexed
        for ii = 1:length(SegData{gg}.Trial) %surgeons
            for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
                for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                    %get the current group, trial, hand, and segment
                    tempstate = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss}.AllStates(:,coluse);
                    
                     %number of points in this iteration
                    [ns,nc] = size(tempstate);
                    
                    if(ii == leaveout)
                        %left out data
                        Xtest = [Xtest; tempstate];
                        Labelstest = [Labelstest; ones(ns,1)*mapping(gg) ];
                    else
                        %training data
                        Xtrain = [Xtrain; tempstate];
                        Labelstrain = [Labelstrain; ones(ns,1)*mapping(gg) ];
                    end
                end
            end
        end
    end
    
    %We have to SubSample to get even counts
    [trainsub,labelsub] = EvenSampleData(Xtrain,Labelstrain);
    
    %train forest
    mdlloo = TreeBagger(trees,trainsub,labelsub,'Method','classification');

    %now we classify with test data
    for gg = grplist

        %grab one classes data
       xtemp = Xtest(Labelstest ==  mapping(gg),:);
       labeltrue = ones(size(xtemp,1),1)*mapping(gg);

       %try classifying 
       esttimecell = predict(mdlloo,xtemp);
       esttime = str2double(esttimecell);
       
       
       %over all
       cesttemp = (mean(esttime) >= 0.5);
       
       
       %cumilative
       cestcumtemp = (CumMean(esttime) >= 0.5);
       
       %score
       scoretimetemp = esttime;

       %find where the last time it changed classification was
       convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
       if(isempty(convtime) )
            convtime = 0;
       end
       
       %store it temporarily
       ClassificationAll = [ClassificationAll; mapping(gg), cesttemp];
       ClassificationTimeAll = [ClassificationTimeAll; labeltrue, cestcumtemp];
       ConvergenceTime = [ConvergenceTime; convtime, mapping(gg)];
       DataStore = [DataStore; xtemp];
       ScoreTime.Class{gg}.run{oo} = scoretimetemp;
    end

end

%convergence time
avgconvtime = mean(ConvergenceTime(:,1))
idxconv1 = find(ConvergenceTime(:,2) == 0);
idxconv2 = find(ConvergenceTime(:,2) == 1);
avgconvtime1 = mean(ConvergenceTime(idxconv1,1))
avgconvtime2 = mean(ConvergenceTime(idxconv2,1))

%inidividual accs (1: liver, 2 : pancreas)
idxTime1 = find(ClassificationTimeAll(:,1) == 0);
idxTime2 = find(ClassificationTimeAll(:,1) == 1);
corrtime1 = ClassificationTimeAll(idxTime1,1) == ClassificationTimeAll(idxTime1,2);
corrtime2 = ClassificationTimeAll(idxTime2,1) == ClassificationTimeAll(idxTime2,2);
acctime1 = mean(corrtime1)
acctime2 = mean(corrtime2)
idx1 = find(ClassificationAll(:,1) == 0);
idx2 = find(ClassificationAll(:,1) == 1);
corr1 = ClassificationAll(idx1,1) == ClassificationAll(idx1,2);
corr2 = ClassificationAll(idx2,1) == ClassificationAll(idx2,2);
acc1 = mean(corr1)
acc2 = mean(corr2)

%overall acc
corrtime = ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2);
acctime = mean(corrtime)
corr = ClassificationAll(:,1) == ClassificationAll(:,2);
acc = mean(corr)

%% Macro accuracy

% macrotime = 82.4278
% macroall = 54.4582

subjects = [29,25,13];
acctime = [81.04,81.4,87.5];
accall = [60.3,50,50];

macrotime = sum(acctime.*subjects) ./ sum(subjects)
macroall = sum(accall.*subjects) ./ sum(subjects)
