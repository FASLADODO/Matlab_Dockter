%% MADE A SMALLER FILE TO JUST LOOK AT INTENT VECTORS

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

%FLS scores
flsScores = cell2mat(DataGlb.content(DataGlb.LogIdx{tsk},flsColumn));

%for access
key.c.PathLength = 1;
key.c.EOM = 2;
key.c.Jerk = 3;
key.c.Crv = 4;
key.c.Acc = 5;
key.c.Time = 6; %not used
key.c.All = {'Path Length','EOM','Jerk','Crv','Acc','Task Time'};

Hands = [1,2];
HandString = {'left', 'right'};


%To avoid having a ton of novice data
allsurgeons = [length(DataGlb.grp.all{tsk}.Idx{myGroups(1)}),length(DataGlb.grp.all{tsk}.Idx{myGroups(2)}),length(DataGlb.grp.all{tsk}.Idx{myGroups(3)})]
min_num_surgeons = min(allsurgeons)
allsurgeons

subset_nov = ArrayIndexSubsets(allsurgeons(1), min_num_surgeons);
subset_int = ArrayIndexSubsets(allsurgeons(2), min_num_surgeons);
subset_exp = ArrayIndexSubsets(allsurgeons(3), min_num_surgeons);

%% Transform each segment to a new coordinate system: Origin is the end of the segment


%for storage
DataAgg = [];
LabelsAgg = [];


%get average segment length
SegmentDuration = [];

%Clear this bad boy
SegData = [];
All = [];

tlength = [];
samplength = [];
totalSegs = [];
totalSurgons = [];

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
    
    
    %randind = sort(randsample(allsurgeons(gg),min_num_surgeons))';
    randind = 1:allsurgeons(gg);
    % Sum up all logs specific group (int, nov, exp)
    for iii = 1:length(randind)
        ii = randind(iii);
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        %Get aggregate metrics for trial
        PathLength_Trial = cell2mat(DataGlb.content(i,toolPathColumn));
        EOM_Trial = cell2mat(DataGlb.content(i,EOMColumn));
        Jerk_Trial = cell2mat(DataGlb.content(i,JerkColumn));
        Crv_Trial = cell2mat(DataGlb.content(i,CrvColumn));
        Acc_Trial = cell2mat(DataGlb.content(i,AccColumn));
        Time_Trial = cell2mat(DataGlb.content(i,timeColumn));
        
        %stash it
        DataAgg{gg}.Trial{iii}(:,key.c.PathLength) = PathLength_Trial ; %to scale for time effects
        DataAgg{gg}.Trial{iii}(:,key.c.EOM) = EOM_Trial;
        DataAgg{gg}.Trial{iii}(:,key.c.Jerk) = Jerk_Trial;
        DataAgg{gg}.Trial{iii}(:,key.c.Crv) = Crv_Trial;
        DataAgg{gg}.Trial{iii}(:,key.c.Acc) = Acc_Trial;
        DataAgg{gg}.Trial{iii}(:,key.c.Time) = Time_Trial;
        
        
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
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.Time = t_seg;
            
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.Pos = [xleft_seg, yleft_seg, zleft_seg ];
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.Vel = [dxleft_seg, dyleft_seg, dzleft_seg ];
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.Acc = [ddxleft_seg, ddyleft_seg, ddzleft_seg ];
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.Jerk = [dddxleft_seg, dddyleft_seg, dddzleft_seg ];
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.Grip = [tempGL ];
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.dGrip = [tempdGL ];
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.VelAngle = [tmpvelmag, tmpvelalpha, tmpvelbeta];
            SegData{gg}.Trial{iii}.Hand{L}.Segment{s}.AccAngle = [tmpaccmag, tmpaccalpha, tmpaccbeta ];
        
            segidx = segidx +1;
            
            SegmentDuration = [SegmentDuration; length(t_seg), gg];
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
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.Time = t_seg;
            
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.Pos = [xright_seg, yright_seg, zright_seg ];
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.Vel = [dxright_seg, dyright_seg, dzright_seg ];
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.Acc = [ddxright_seg, ddyright_seg, ddzright_seg ];
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.Jerk = [dddxright_seg, dddyright_seg, dddzright_seg ];
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.Grip = [tempGR ];
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.dGrip = [tempdGR ];
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.VelAngle = [tmpvelmag, tmpvelalpha, tmpvelbeta];
            SegData{gg}.Trial{iii}.Hand{R}.Segment{s}.AccAngle = [tmpaccmag, tmpaccalpha, tmpaccbeta ];
        
            segidx = segidx +1;
            
            SegmentDuration = [SegmentDuration; length(t_seg), gg];
        end
        
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
    totalSurgons(gg) = iii;
end
fprintf('Done \n');

totalSegs
totalSurgons


durationNovice = SegmentDuration(SegmentDuration(:,2) == 1,1);
durationExpert = SegmentDuration(SegmentDuration(:,2) == 3,1);

mean(durationNovice)
std(durationNovice)
mean(durationExpert)
std(durationExpert)

%% Lets look at INTENT VECTORS leave one out
% Motion efficiency made you think of this
% How often does the person change direction from their overall objective


%For suturing we get
% Novice: 100% classification, Expert: 96%
%For cuttin

%which states we want to keep
StrLabels = {'Novice';'Intermediate';'Expert'};
grplist = [1,3];
ploton = false;
ThreshStore = [];
ClassificationEstimate = [];
ClassificationTrue = [];
stashSurgeon =[];
allSegmentCounts = [];

saveVars = [];

%Store counts for later
CountsStore = [];

%All possible keepers
idxall = 1:min_num_surgeons;
ModelUse = 3; %experts

fprintf('Computing Intent Vectors leave one out... \n');

%sublim = size(subset_nov,1);
%sublim = max(allsurgeons);
sublim = allsurgeons(1);
for ff = 1:1
    
%     shiftidx(1,:) = subset_nov(ff,:);
%     %shiftidx(2,:) = subset_int(ff,:);
%     shiftidx(3,:) = idxall;
%     shiftidx
    
    %Loop through which one subject will be left out
    for oo = 1:sublim
    %for oo = idxall
        %fprintf('shift: %d of %d\n',ff,sublim);
        fprintf('epoch: %d of %d\n',oo,sublim);

        %choose what we're leaving out
%         idxtrain = idxall;
%         idxtrain(oo) = [];
%         idxtrain 
%         shiftidx(1,idxtrain)
%         idxtest = idxall(oo);
        for glo = grplist
            idxtrain{glo} = 1:allsurgeons(glo);
            glo
            out = mod(oo,allsurgeons(glo));
            if(out == 0)
                out = allsurgeons(glo);
            end
            idxtrain{glo}(out) = [];
            idxtest{glo} = out;
        end

        %clear this monster
        SegTrain = [];
        SegTest = [];

        DataTrain = [];
        DataTest = [];
        LabelsTrain = [];
        LabelsTest = [];

        %aggregate metrics
        AggTrain = [];
        AggTest = [];
        
        % loop through samples to compute intent vector info
        for gg = grplist %skill levels
            %for smapling training vs validation
            surgeonid = 1;
            idtrain = 1;
            idtest = 1;

            
            for ii = 1:length(SegData{gg}.Trial) %surgeons
            %for iii = idxall %surgeons
                %ii = shiftidx(gg,iii);
                
                segid = 1;
                for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
                    for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                        %get the current group, trial, hand, and segment
                        struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};

                        %XYZ data
                        Data = struct.Pos;
                        [NN,SS] = size(Data);

                        if(NN < 5)
                            %skip it
                            continue;
                        end
                        %get intent vector angles and progress
                        [IV_Angle,endp] = IntentVectorSegments(Data);
                        [Progress] = IntentVectorProgress(Data);

                        %stash it
                        dtemp = [IV_Angle, Progress];

                        if(any(surgeonid==idxtrain{gg}))
                            DataTrain = [DataTrain; dtemp];
                            LabelsTrain = [LabelsTrain; ones(NN,1)*gg];
                            SegTrain{gg}.Trial{idtrain}.Segment{segid} = dtemp;
                            %aggregate measures
                            AggTrain{gg}.Trial{idtrain} = DataAgg{gg}.Trial{ii};
                        elseif(any(surgeonid==idxtest{gg}))
                            DataTest = [DataTest; dtemp];
                            LabelsTest = [LabelsTest; ones(NN,1)*gg];
                            SegTest{gg}.Trial{idtest}.Segment{segid} = dtemp;
                            %aggregate measures
                            AggTest{gg}.Trial{idtest} = DataAgg{gg}.Trial{ii};
                        end

                        segid = segid + 1;
                    end
                end
                if(any(surgeonid==idxtrain{gg}))
                    idtrain = idtrain + 1;
                elseif(any(surgeonid==idxtest{gg}))
                    idtest = idtest + 1;
                end
                surgeonid = surgeonid + 1;
            end

        end

        %Now get seperability and gaussian model from 'good data'
        [Difference,ClassData,ProbData] = SimpleRelativeRBFTrain(DataTrain,LabelsTrain);
        [Model,GoodData,threshG] = SimpleRelativeRBFProbs(Difference,ClassData,ProbData);


        %now train the threshold for within expert probability or not (T_P)
        stank = SegTrain;
        
        %for combined
        CombineFeatures = [];
        CombineLabels = [];

        expmodel = Model{ModelUse};
        stashprobs = [];
        stashSurgeonT =[];
        for gg = grplist

             for ii = 1:length(stank{gg}.Trial) %surgeons
                 perSurgeonCount = [];
                 for ss = 1:length(stank{gg}.Trial{ii}.Segment) %segments
                     dtemp = stank{gg}.Trial{ii}.Segment{ss};

                     %get all IV probabilities for this segment
                     probtemp = mvnpdf(dtemp,expmodel.mean,expmodel.sigma); 

                     %check how many parts of segment are below expert thresh
                     class = probtemp < expmodel.thresh;

                     %stash that
                     stashprobs = [ stashprobs; sum(class), gg];

                     perSurgeonCount = [perSurgeonCount; sum(class)];
                     
                 end

                 CombineFeatures = [CombineFeatures; mean(perSurgeonCount), AggTrain{gg}.Trial{ii}];
                 CombineLabels = [CombineLabels; gg ];
             end
        end

        
        %Collect all features in a qda
        MdlLinear = fitcdiscr(CombineFeatures,CombineLabels);
        

        % now try classifying test data with the expert probability model
        stank = SegTest;

        expmodel = Model{ModelUse};
        stashprobs = [];
        for gg = grplist

             for ii = 1:length(stank{gg}.Trial) %surgeons
                 perSurgeonCount = [];
                 for ss = 1:length(stank{gg}.Trial{ii}.Segment)
                     dtemp = stank{gg}.Trial{ii}.Segment{ss};

                     %get all IV probabilities for this segment
                     probtemp = mvnpdf(dtemp,expmodel.mean,expmodel.sigma); 

                     %check how many parts of segment are below expert thresh
                     class = probtemp < expmodel.thresh;

                     %stash that
                     stashprobs = [ stashprobs; sum(class), gg];

                     perSurgeonCount = [perSurgeonCount; sum(class)];
                 end
                 
                 %leaver = shiftidx(gg,idxtest)
                 leaver = idxtest{gg}
                 
                 %get combine vector for this guy
                 combtemp = [mean(perSurgeonCount), AggTest{gg}.Trial{ii}];
                 
                 %test classify
                 estLabel = predict(MdlLinear,combtemp);
                 
                 %store the counts for this particular surgeon
                 stashSurgeon{gg}(leaver) = estLabel;

                 %check how we classify the leave one out given the number of
                 ClassificationEstimate = [ClassificationEstimate; estLabel];
                 ClassificationTrue = [ClassificationTrue; gg];
                 
             end

        end
       

    end
    
end
fprintf('Done \n');

ClassificationEstimate
corr = ClassificationEstimate == ClassificationTrue; 
acctotal = mean(corr)

%combined accs:
%(confirmed)
%peg transfer (1): 93.1% total, 97.6% macro 100% novices, 86.2% Expert 29 nov, 6 exp
%(confirmed)
%pattern cutting: (2) 98% total, 97.2% macro, 96% nov, 100% exp 25 nov, 10 experts
%(confirmed)
%suturing (3): 96.2% total, 95.2% macro, 92.3% nov, 100% exp, 13 nov, 8 Experts


%only do each surgeon once
maxnov = allsurgeons(1);
maxexp = allsurgeons(3);

%stash only nov and exp separately
novStore = ClassificationEstimate(ClassificationTrue == 1,:);
expStore = ClassificationEstimate(ClassificationTrue == 3,:);
novSub = novStore(1:maxnov,:);
expSub = expStore(1:maxexp,:);
novcorr = novSub == 1;
expcorr = expSub == 3;
novacc = mean(novcorr)
expacc = mean(expcorr)

%macro accuracy
accz = [novacc,expacc];
countz = [allsurgeons(1),allsurgeons(3)];
macroacc = sum( (accz.*countz) ) ./ sum(countz)


return;

%% get each classes estimates

%stash only nov and exp separately
novStore = CountsStore(CountsStore(:,3) == 1,:);
expStore = CountsStore(CountsStore(:,3) == 3,:);

%only do each surgeon once
maxnov = max(novStore(:,5))
maxexp = max(expStore(:,5))

%only first loop all indices
novSub = novStore(1:maxnov,:);
expSub = expStore(1:maxexp,:);
allSub = [novSub; expSub];

%get corrects 
novCorr = novSub(:,3) == novSub(:,4);
expCorr = expSub(:,3) == expSub(:,4);
corrAll = allSub(:,3) == allSub(:,4);

%get accuracy for each
novAcc = mean(novCorr)
expAcc = mean(expCorr)
accreal = mean(corrAll)

%% get total segments

novEst = CountsStore(CountsStore(:,4) == 1,:);
expEst = CountsStore(CountsStore(:,4) == 3,:);
fsize = 14;
plotstr = {'Novice';'Intermediate';'Expert'};

figure
hax=axes; 
h1=scatter(novEst(:,3), novEst(:,1),60, 'ro')
hold on
h2=scatter(expEst(:,3), expEst(:,1),60, 'co')
for ppp = 1:length(ThreshStore)
    hold on
    h3=line(get(hax,'XLim'),[ThreshStore(ppp) ThreshStore(ppp)],'Color',[0 0 0])
end
hold off
hs = legend([h1(1),h2(1),h3(1)],'Est. Novice','Est. Expert','Threshold')
str = sprintf(' C_{g} Total Counts (Task: %s), (Accuracy %.2f)',TaskStr{tsk},accreal);
%title(str)
xlabel('True Skill Level','FontSize',fsize)
ylabel('Total Demerits (SK)','FontSize',fsize)
xlim([0.8, 3.2])
set(gca,'XTick',[1,3],'XTickLabel',{'Novice';'Expert'})
hs.FontSize = fsize;
hs.Location = 'NorthEast';


%% plot average segment counts

allSegNov = allSegmentCounts(allSegmentCounts(:,1) == 1,2);
allSegExp = allSegmentCounts(allSegmentCounts(:,1) == 3,2);

mean(allSegNov)
std(allSegNov)


mean(allSegExp)
std(allSegExp)

% ans =
%    65.9194
% ans =
%   105.2392
% ans =
%    22.5568
% ans =
%    27.6695


plotstr = {'Novice';'Intermediate';'Expert'};
fsize = 14;
figure
boxplot(allSegmentCounts(:,2),plotstr(allSegmentCounts(:,1)),'Notch','on')
xlabel('Skill Level','FontSize',fsize);
ylabel('Segment Demerits (y_{i})','FontSize',fsize);
%title('Segment Demerits (y_{i}) ','FontSize',fsize);
set(gca,'YScale','log')

%% Figure out total overall accuracy

accz = [96.5,83,100,100,100,92.3]
countz = [29,6,25,10,13,8]

overallacc = sum( (accz.*countz) ./ sum(countz))

accpegtx = sum( (accz(1:2).*countz(1:2)) ./ sum(countz(1:2)))
acccut = sum( (accz(3:4).*countz(3:4)) ./ sum(countz(3:4)))
accsut = sum( (accz(5:6).*countz(5:6)) ./ sum(countz(5:6)))

%%

return


