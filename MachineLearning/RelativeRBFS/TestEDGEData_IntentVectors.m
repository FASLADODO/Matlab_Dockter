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
    
    
    randind = sort(randsample(allsurgeons(gg),min_num_surgeons))'
    % Sum up all logs specific group (int, nov, exp)
    for iii = 1:allsurgeons(gg) %length(randind)
        %ii = randind(iii);
        ii = iii;
        
        i = DataGlb.grp.all{tsk}.Idx{myGroups(gg)}(ii);
        
        %Get FLS score for that
        FLS_Trial = cell2mat(DataGlb.content(i,flsColumn));
        SegData{gg}.Trial{iii}.FLSScore = FLS_Trial;
        
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




%% INTENT VECTORS SAMPLE PLOTS

StrLabels = {'Novice','Expert'};
grplist = [1,3];
fsize = 16;

IDX_NOV = 103; %240; %randi(totalSegs(1),1) %103 or 101?
IDX_EXP = -1; %30; %30 %randi(totalSegs(3),1)
% loop through samples
for gg = grplist %skill levels
    segid = 1;
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};
                
                Data = struct.Pos;
                [IV_Angle,endp,IV_Diff] = IntentVectorSegments(Data);
                [Progress] = IntentVectorProgress(Data);
                
                if(segid == IDX_NOV && gg == 1)
                    figure(1)
                    scatter3(Data(:,1),Data(:,2),Data(:,3),40,IV_Angle,'filled')
                    hold on
                    plot3(endp(:,1),endp(:,2),endp(:,3),'r','LineWidth',1);
                    hold off
                    xlabel('x (cm)','FontSize',fsize)
                    ylabel('y (cm)','FontSize',fsize)
                    zlabel('z (cm)','FontSize',fsize)
                    %title('Intent Vector Angle','FontSize',fsize)
                    %title('IV ANGLE NOV')
                    colormap winter
                    hc = colorbar;
                    ylabel(hc, 'IVA','FontSize',fsize)
                    
                    figure(2)
                    scatter3(Data(:,1),Data(:,2),Data(:,3),40,Progress,'filled')
                    hold on
                    plot3(endp(:,1),endp(:,2),endp(:,3),'r','LineWidth',1);
                    hold off
                    xlabel('x (cm)','FontSize',fsize)
                    ylabel('y (cm)','FontSize',fsize)
                    zlabel('z (cm)','FontSize',fsize)
                    %title('Intent Vector Progress','FontSize',fsize)
                    %title('IV PROGRESS NOV')
                    colormap winter
                    hc = colorbar;
                    ylabel(hc, 'IVP','FontSize',fsize)
                    
%                     figure(5)
%                     scatter3(Data(:,1),Data(:,2),Data(:,3),40,IV_Diff)
%                     hold on
%                     plot3(endp(:,1),endp(:,2),endp(:,3),'g');
%                     hold off
%                     xlabel('x (cm)','FontSize',fsize)
%                     ylabel('y (cm)','FontSize',fsize)
%                     zlabel('z (cm)','FontSize',fsize)
%                     title('IV ANGLE DIFF')
%                     colormap cool
%                     hc = colorbar;
%                     ylabel(hc, 'dIVA','FontSize',fsize)
                end
                if(segid == IDX_EXP && gg == 3)
                    figure(3)
                    scatter3(Data(:,1),Data(:,2),Data(:,3),40,IV_Angle,'filled')
                    hold on
                    plot3(endp(:,1),endp(:,2),endp(:,3),'r','LineWidth',1);
                    hold off
                    xlabel('x (cm)','FontSize',fsize)
                    ylabel('y (cm)','FontSize',fsize)
                    zlabel('z (cm)','FontSize',fsize)
                    %title('Intent Vector Angle','FontSize',fsize)
                    %title('IV ANGLE EXP','FontSize',fsize)
                    colormap winter
                    hc = colorbar;
                    ylabel(hc, 'IVA','FontSize',fsize)
                    
                    figure(4)
                    scatter3(Data(:,1),Data(:,2),Data(:,3),40,Progress,'filled')
                    hold on
                    plot3(endp(:,1),endp(:,2),endp(:,3),'r','LineWidth',1);
                    hold off
                    xlabel('x (cm)','FontSize',fsize)
                    ylabel('y (cm)','FontSize',fsize)
                    zlabel('z (cm)','FontSize',fsize)
                    %title('Intent Vector Progress','FontSize',fsize)
                    %title('IV PROGRESS EXP')
                    colormap winter
                    hc = colorbar;
                    ylabel(hc, 'IVP','FontSize',fsize)
                    
%                     figure(6)
%                     scatter3(Data(:,1),Data(:,2),Data(:,3),40,IV_Diff)
%                     hold on
%                     plot3(endp(:,1),endp(:,2),endp(:,3),'g');
%                     hold off
%                     xlabel('x (cm)','FontSize',fsize)
%                     ylabel('y (cm)','FontSize',fsize)
%                     zlabel('z (cm)','FontSize',fsize)
%                     title('IV ANGLE DIFF')
%                     colormap cool
%                     hc = colorbar;
%                     ylabel(hc, 'dIVA','FontSize',fsize)
                end
                
                segid = segid + 1;
            end
        end
    end;
end

%% Lets look at INTENT VECTORS
% Motion efficiency made you think of this
% How often does the person change direction from their overall objective


%which states we want to keep
StrLabels = {'Novice','Expert'};
grplist = [1,3];

%clear this monster
SegTrain = [];
SegValidate = [];
SegTest = [];

DataTrain = [];
DataValidate = [];
DataTest = [];
LabelsTrain = [];
LabelsValidate = [];
LabelsTest = [];


%Manual (dont try this at home kids)
% idxtrain = [4,5,6,7,8];
% idxval = [1,2];
% idxtest = [3];
idxtrain = [1,2,3,4,5,6,7,8];
idxval = [];
idxtest = [];

%sample into training/validation/testing
trainRatio = 0.625;
valRatio = 0.25;
testRatio = 0.125;

fprintf('Computing Intent Vectors...');
% loop through samples
for gg = grplist %skill levels
    %for smapling training vs validation
    %[idxtrain, idxval, idxtest] = dividerand(totalSurgons(gg),trainRatio,valRatio,testRatio)
    surgeonid = 1;
    idtrain = 1;
    idval = 1;
    idtest = 1;
    
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        
        segid = 1;
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};
                
                %XYZ data
                Data = struct.Pos;
                [NN,SS] = size(Data);
                
                if(NN < 4)
                   continue; 
                end

                %get intent vector angles and progress
                [IV_Angle,endp,IV_Diff] = IntentVectorSegments(Data);
                [Progress] = IntentVectorProgress(Data);
                
                %stash it
                dtemp = [IV_Angle, Progress];
                %dtemp = [IV_Angle, Progress, IV_Diff];
                
                if(any(surgeonid==idxtrain))
                    DataTrain = [DataTrain; dtemp];
                    LabelsTrain = [LabelsTrain; ones(NN,1)*gg];
                    SegTrain{gg}.Trial{idtrain}.Segment{segid} = dtemp;
                elseif(any(surgeonid==idxval))
                    DataValidate = [DataValidate; dtemp];
                    LabelsValidate = [LabelsValidate; ones(NN,1)*gg];
                    SegValidate{gg}.Trial{idval}.Segment{segid} = dtemp;
                elseif(any(surgeonid==idxtest))
                    DataTest = [DataTest; dtemp];
                    LabelsTest = [LabelsTest; ones(NN,1)*gg];
                    SegTest{gg}.Trial{idtest}.Segment{segid} = dtemp;
                end
                
                segid = segid + 1;
            end
        end
        if(any(surgeonid==idxtrain))
            idtrain = idtrain + 1;
        elseif(any(surgeonid==idxval))
            idval = idval + 1;
        elseif(any(surgeonid==idxtest))
            idtest = idtest + 1;
        end
        surgeonid = surgeonid + 1;
    end
    
end
fprintf('Done \n');


labstr = {'Novice';'Intermediate';'Expert'}
figure
%gscatter(DataTrain(:,1),DataTrain(:,2),labstr(LabelsTrain));
%gscatter3(DataTrain(:,1),DataTrain(:,2),DataTrain(:,3),LabelsTrain);
hs=gscatter(DataTrain(:,1),DataTrain(:,2),labstr(LabelsTrain),'rc');
hold off
xlabel('IVA (rad)','FontSize',fsize)
ylabel('IVP ','FontSize',fsize)
% xlabel('IVA (rad)','FontSize',fsize)
% ylabel('dIVA ','FontSize',fsize)
%title('Intent Vectors','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';
xlim([0,pi])
ylim([-2.2,3.2])

%check our ratios
allsizes = [length(SegTrain{gg}.Trial),length(SegValidate{gg}.Trial),length(SegTest{gg}.Trial)];
sumsizes = sum(allsizes);
actualratios = allsizes ./ sumsizes


%% FOR PLOTTING TRUE EXPERT AND W_EXP (Fig 4 and 5) (PAPER!)

[Difference,ClassData,Model] = SimpleRelativeRBFTrain(DataTrain,LabelsTrain);
[ProbModel,GoodData,threshG] = SimpleRelativeRBFProbs(Difference,ClassData,Model);

cc = 3;
allProbs = mvnpdf(ClassData{cc},ProbModel{cc}.mean,ProbModel{cc}.sigma);

figure
% gscatter(DataTrain(:,1),DataTrain(:,2),LabelsTrain);
cc = 1;
h1 = scatter(ClassData{cc}(:,1),ClassData{cc}(:,2),'r.')
hold on
cc = 3;
h2 = scatter(ClassData{cc}(:,1),ClassData{cc}(:,2),'c.')
hold on
cc = 3;
h3 = scatter(Model{cc}.gooddata(:,1),Model{cc}.gooddata(:,2),'bo')
% h3 = scatter(GoodData{cc}(:,1),GoodData{cc}(:,2),'bo')
hold off
legend([h1(1),h2(1),h3(1)],'Novice','Expert','True Expert')
xlabel('IVA (rad)','FontSize',fsize)
ylabel('IVP','FontSize',fsize)
%title('True Expert Data','FontSize',fsize)
lgz= findobj(gcf,'tag','legend'); 
set(lgz,'FontSize',fsize)


figure
cc = 1;
h1 = scatter(ClassData{cc}(:,1),ClassData{cc}(:,2),'r.')
hold on
cc = 3;
h2 = scatter3(ClassData{cc}(:,1),ClassData{cc}(:,2),ones(length(ClassData{cc}),1)*0.001,'c.')
hold on
cc = 3;
h3 = Surface3D(Model{cc}.gooddata(:,1),Model{cc}.gooddata(:,2),Model{cc}.prob,'mesh',[0,pi;-2.2,3.2])
hold off
legend([h1(1),h2(1)],'Novice','Expert')
xlabel('IVA (rad)','FontSize',fsize)
ylabel('IVP','FontSize',fsize)
zlabel('W_{exp}','FontSize',fsize)
%title('Relevance Weights','FontSize',fsize)
hc = colorbar;
ylabel(hc, 'W_{exp}','FontSize',fsize)
lgz= findobj(gcf,'tag','legend'); 
set(lgz,'FontSize',fsize)
xlim([0,pi])
ylim([-2.2,3.2])


figure
cc = 1;
h1 = scatter(ClassData{cc}(:,1),ClassData{cc}(:,2),'r.')
hold on
cc = 3;
h2 = scatter(ClassData{cc}(:,1),ClassData{cc}(:,2),'c.')
hold on
cc = 3;
h3 = Surface3D(DataTrain(:,1),DataTrain(:,2),allProbs)
hold off
legend([h1(1),h2(1)],'Novice','Expert')
xlabel('IVA (rad)','FontSize',fsize)
ylabel('IVP','FontSize',fsize)
zlabel('W_{exp}','FontSize',fsize)
%title('Relevance Weights','FontSize',fsize)
hc = colorbar;
ylabel(hc, 'W_{exp}','FontSize',fsize)
lgz= findobj(gcf,'tag','legend'); 
set(lgz,'FontSize',fsize)
xlim([0,pi])
ylim([-2.2,3.2])


%% determine RBF for single class

CUSE = 3;
DifferenceTrain = RelativeRBFSingleClass(DataTrain,LabelsTrain,CUSE);

figure
gscatter(DataTrain(:,1),DataTrain(:,2),labstr(LabelsTrain))
hold on
Surface3D(DataTrain(:,1),DataTrain(:,2),DifferenceTrain);
hold off
title('rbf class exp')
%zlabel('$W_{rbf}$','FontSize',fsize,'interpreter','latex')
hc = colorbar;
ylabel(hc, 'W_{rbf}','FontSize',fsize)
view ([40,30])
%axis([-3 8 -3 8 0 1.2])

%% Compute optimal threshold from information gain

Direction = [3;1];
conservative = 2; 

[ThreshTrain,IG,Stash] = selectThresholdIG(DifferenceTrain,LabelsTrain,Direction);
ThreshTrain = ThreshTrain*conservative
IG
threshold = 0.1*prctile(DifferenceTrain,99)

figure
scatter(Stash(:,1),Stash(:,2))
xlabel('threshold')
ylabel('IG')
title('threshold choice from IG')

%%  estimate class and plot threshdold

[classest,sep] = RelativeRBFSingleClassOnline(DataTrain,LabelsTrain,DataTrain,Direction,threshold);


figure
gscatter(DataTrain(:,1),DataTrain(:,2),classest)
title('class estimates')

zh = threshold;
mu = [1.5,1];
sq = [2,3];
figure
gscatter(DataTrain(:,1),DataTrain(:,2),classest)
hold on
plotZplane(mu,sq,zh);
hold on
Surface3D(DataTrain(:,1),DataTrain(:,2),DifferenceTrain);
hold off
title('IG Threshold')

corr = classest == LabelsTrain;

acc = mean(corr)


%% train counts threshold

stank = SegTrain;

randseg = 50;


stashprobs = [];
stashSurgeon =[];
for gg = grplist
     segid = 1;
     
     for ii = 1:length(stank{gg}.Trial) %surgeons
         perSurgeonCount = [];
         for ss = 1:length(stank{gg}.Trial{ii}.Segment)
             dtemp = stank{gg}.Trial{ii}.Segment{ss};
             
             %get all IV probabilities for this segment
             [class,sep] = RelativeRBFSingleClassOnline(DataTrain,LabelsTrain,dtemp,Direction,threshold);
             
             %stash that
             stashprobs = [ stashprobs; sum(class), gg];
             
             perSurgeonCount = [perSurgeonCount; sum(class)];
             
             segid = segid + 1;
             if(segid == randseg && gg == 1)
                figure
                gscatter(DataTrain(:,1),DataTrain(:,2),classest)
                hold on
                Surface3D(dtemp(:,1),dtemp(:,2),sep);
                hold off
                title('IG NOV')
             elseif(segid == randseg && gg == 3)
                figure
                gscatter(DataTrain(:,1),DataTrain(:,2),classest)
                hold on
                Surface3D(dtemp(:,1),dtemp(:,2),sep);
                hold off
                title('IG EXP')

             end
         end
         
         stashSurgeon{gg}(ii) = mean(perSurgeonCount);
     end
     
end

figure
gscatter(stashprobs(:,2), stashprobs(:,1), stashprobs(:,2))
title('all segments counts')

figure
cc = 1;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'ro')
hold on
cc = 2;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'go')
hold on
cc = 3;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'bo')
hold off
title('total for surgeons')


%get thresh
XCountData = [stashSurgeon{1}'; stashSurgeon{3}'];
LabelXCount = [ones(length(stashSurgeon{1}),1)*1;ones(length(stashSurgeon{3}),1)*3];
WCount = ones(length(XCountData),1)/length(XCountData);
[ThreshC,score] = DecisionStumpBasic(XCountData,WCount,LabelXCount);
ThreshC

%% now try classifying each segment with model

stank = SegTest;

stashprobs = [];
stashSurgeon =[];
for gg = grplist
    
     for ii = 1:length(stank{gg}.Trial) %surgeons
         perSurgeonCount = [];
         for ss = 1:length(stank{gg}.Trial{ii}.Segment)
             dtemp = stank{gg}.Trial{ii}.Segment{ss};
             
             %get all IV probabilities for this segment
             [class,sep] = RelativeRBFSingleClassOnline(DataTrain,LabelsTrain,dtemp,Direction,threshold);
             
             %stash that
             stashprobs = [ stashprobs; sum(class), gg];
             
             perSurgeonCount = [perSurgeonCount; sum(class)];
         end
         
         stashSurgeon{gg}(ii) = mean(perSurgeonCount);
     end
     
end

figure
gscatter(stashprobs(:,2), stashprobs(:,1), stashprobs(:,2))
title('all segments counts')

figure
cc = 1;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'ro')
cc = 3;
hold on
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'bo')
hold off
title('total for surgeons')






%%
return




%% train threshold

stank = SegTrain;

expmodel = Model{3};
stashprobs = [];
stashSurgeon =[];
for gg = grplist
    
     for ii = 1:length(stank{gg}.Trial) %surgeons
         perSurgeonCount = [];
         for ss = 1:length(stank{gg}.Trial{ii}.Segment)
             dtemp = stank{gg}.Trial{ii}.Segment{ss};
             
             %get all IV probabilities for this segment
             probtemp = RBFClassProbability(expmodel.gooddata,dtemp);
             %probtemp = mvnpdf(dtemp,expmodel.mean,expmodel.sigma); 
             
             %check how many parts of segment are below expert thresh
             class = probtemp < expmodel.thresh;
             
             %stash that
             stashprobs = [ stashprobs; sum(class), gg];
             
             perSurgeonCount = [perSurgeonCount; sum(class)];
         end
         
         stashSurgeon{gg}(ii) = mean(perSurgeonCount);
     end
     
end

figure
gscatter(stashprobs(:,2), stashprobs(:,1), stashprobs(:,2))
title('all segments counts')

figure
cc = 1;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'ro')
hold on
cc = 2;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'go')
hold on
cc = 3;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'bo')
hold off
title('total for surgeons')


%get thresh
XCountData = [stashSurgeon{1}'; stashSurgeon{3}'];
LabelXCount = [ones(length(stashSurgeon{1}),1)*1;ones(length(stashSurgeon{3}),1)*3];
WCount = ones(length(XCountData),1)/length(XCountData);
[ThreshC,score] = DecisionStumpBasic(XCountData,WCount,LabelXCount);
ThreshC

%% now try classifying each segment with model

stank = SegTest;

expmodel = Model{2};
stashprobs = [];
stashSurgeon =[];
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
         
         stashSurgeon{gg}(ii) = mean(perSurgeonCount);
     end
     
end

figure
gscatter(stashprobs(:,2), stashprobs(:,1), stashprobs(:,2))
title('all segments counts')

figure
cc = 1;
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'ro')
cc = 3;
hold on
scatter(ones(length(stashSurgeon{cc}),1)*cc,stashSurgeon{cc},'bo')
hold off
title('total for surgeons')


%%

return






%% Old Code


Data_IV = [];
Data_Changes = [];
Labels_IV = [];
Labels_Changes = [];

maxsegsnov = totalSegs(1);
maxsegsexp= totalSegs(3);
IDX_NOV = randi(maxsegsnov,1)
IDX_EXP = randi(maxsegsexp,1)
% loop through samples
for gg = grplist %skill levels
    seglength = [];
    segid = 1;
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};
                
                %XYZ data
                Data = struct.Pos;
                [NN,SS] = size(Data);

                %get intent vector angles and progress
                [IV_Angle,endp] = IntentVectorSegments(Data);
                [Progress] = IntentVectorProgress(Data);
                

                %now look at how many times we change direction
                IV_Signed = pi/2 - IV_Angle;
                IV_Changes = abs(diff(sign(IV_Signed))) > eps;
                IV_Changes(end+1) = IV_Changes(end);
                Count_Changes = sum(IV_Changes);

                %stash it
                %Data_IV = [Data_IV; IV_Angle, Progress, ones(NN,1)*Count_Changes];
                Data_IV = [Data_IV; IV_Angle, Progress];
                Labels_IV = [ Labels_IV; ones(NN,1)*gg];
                Data_Changes = [Data_Changes; Count_Changes];
                Labels_Changes = [ Labels_Changes; gg];
                
                if(segid == IDX_NOV && gg == 1)
                    figure
                    scatter3(Data(:,1),Data(:,2),Data(:,3),20,Progress)
                    hold on
                    plot3(endp(:,1),endp(:,2),endp(:,3),'g');
                    hold off
                    str = sprintf('IV nov, changes = %d', Count_Changes);
                    title(str)
                    colormap cool
                    colorbar
                end
                if(segid == IDX_EXP && gg == 3)
                    figure
                    scatter3(Data(:,1),Data(:,2),Data(:,3),20,Progress)
                    hold on
                    plot3(endp(:,1),endp(:,2),endp(:,3),'g');
                    hold off
                    str = sprintf('IV exp, changes = %d', Count_Changes);
                    title(str)
                    colormap cool
                    colorbar
                end
                
                seglength(segid) = size(Data,1);
                segid = segid + 1;
            end
        end
    end
    %segid
end


%% For plotting

%look at average changes
CHANGES_NOV = Data_Changes(Labels_Changes == 1,1);
CHANGES_EXP = Data_Changes(Labels_Changes == 3,1);
CHANGESNOV_AVG = mean(CHANGES_NOV)
CHANGESEXP_AVG = mean(CHANGES_EXP)

%look at average IV angle
IV_NOV = Data_IV(Labels_IV == 1,1);
IV_EXP = Data_IV(Labels_IV == 3,1);
IVNOV_AVG = mean(IV_NOV);
IVEXP_AVG = mean(IV_EXP);

%percentage over 180 backwards
ERROR_NOV = IV_NOV > pi/2;
ERROR_EXP = IV_EXP > pi/2;

sum(ERROR_NOV)
sum(ERROR_EXP)

if(true)
    figure
    gscatter(Data_IV(:,1),Data_IV(:,2),Labels_IV)
    hold on
    %scatter(1,CHANGESNOV_AVG,'go')
    hold on
    %scatter(3,CHANGESEXP_AVG,'go')
    hold off
    xlabel('IV Angle')
    ylabel('motion progress')
end

[Difference,ClassData,ProbData] = SimpleRelativeRBFTrain(Data_IV,Labels_IV);
[Model,GoodData] = SimpleRelativeRBFProbs(Difference,ClassData,ProbData);

figure
gscatter(Data_IV(:,1),Data_IV(:,2),Labels_IV);
hold on
cc = 2;
scatter(GoodData{cc}(:,1),GoodData{cc}(:,2),'bo')
hold off
title('good data exp')

figure
gscatter(Data_IV(:,1),Data_IV(:,2),Labels_IV);
hold on
cc = 2;
scatter(GoodData{cc}(:,1),GoodData{cc}(:,2),'bo')
hold on
Surface3D(GoodData{cc}(:,1),GoodData{cc}(:,2),Model{cc}.prob);
hold off
title('good data exp w/ probz')

figure
gscatter(Data_IV(:,1),Data_IV(:,2),Labels_IV);
hold on
cc = 1;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),ProbData{cc}(:,1));
hold off
title('probs class 1')

figure
gscatter(Data_IV(:,1),Data_IV(:,2),Labels_IV);
hold on
cc = 2;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),ProbData{cc}(:,1));
hold off
title('probs class 3')



%% Lets look at motion efficiency


Data_M = [];

maxsegsnov = totalSegs(1);
maxsegsexp= totalSegs(3);
IDX_NOV = randi(maxsegsnov,1)
IDX_EXP = randi(maxsegsexp,1)
% loop through samples
for gg = grplist %skill levels
    seglength = [];
    segid = 1;
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};
                
                Data = struct.Pos;
                [ME,X1,X2,params] = MotionEfficieny(Data);
                
                Data_M = [Data_M; ME, ones(length(ME),1)*gg];
                
                if(segid == IDX_NOV && gg == 1)
                    figure
                    scatter3(Data(:,1),Data(:,2),Data(:,3),20,ME)
                    hold on
                    endp = [X1;X2];
                    plot3(endp(:,1),endp(:,2),endp(:,3),'g');
                    hold off
                    title('motion efficiency nov')
                    colormap cool
                    colorbar
                end
                if(segid == IDX_EXP && gg == 3)
                    figure
                    scatter3(Data(:,1),Data(:,2),Data(:,3),20,ME)
                    hold on
                    endp = [X1;X2];
                    plot3(endp(:,1),endp(:,2),endp(:,3),'g');
                    hold off
                    title('motion efficiency exp')
                    colormap cool
                    colorbar
                end
                
                seglength(segid) = size(Data,1);
                segid = segid + 1;
            end
        end
    end
    [M,I] = max(seglength);
end

%look at average motion efficiencies
ME_NOV = Data_M(Data_M(:,2) == 1,1);
ME_EXP = Data_M(Data_M(:,2) == 3,1);
MENOV_AVG = mean(ME_NOV);
MEEXP_AVG = mean(ME_EXP);

if(false)
    figure
    gscatter(Data_M(:,2),Data_M(:,1),Data_M(:,2))
    hold on
    scatter(1,MENOV_AVG,'go')
    hold on
    scatter(3,MEEXP_AVG,'go')
    hold off
    xlabel('class')
    ylabel('efficiency')
end


%% Lets look at autocorrelation? not great

kac = 6;

Data_AC = [];

% loop through samples
for gg = grplist %skill levels
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};
                
                PosTemp = struct.VelAngle(:,1);
                VelTemp = struct.Acc(:,3);
                [Rpos,gp] = autocorrelation_simple(PosTemp,kac);
                [Rvel,gv] = autocorrelation_simple(VelTemp,kac);
                if(gp && gv )
                    Data_AC = [Data_AC; Rpos, Rvel, ones(length(Rpos),1)*gg];
                end
            end
        end
    end
end

figure
gscatter(Data_AC(:,1),Data_AC(:,2),Data_AC(:,3))
xlabel('posac')
ylabel('velac')


%% Lets look at FFT? 

Data_AC = [];

% loop through samples
for gg = grplist %skill levels
    for ii = 1:length(SegData{gg}.Trial) %surgeons
        for hh = 1:length(SegData{gg}.Trial{ii}.Hand) %hands
            for ss = 1:length(SegData{gg}.Trial{ii}.Hand{hh}.Segment)
                %get the current group, trial, hand, and segment
                struct = SegData{gg}.Trial{ii}.Hand{hh}.Segment{ss};
                
                PosTemp = struct.VelAngle(:,1);
                VelTemp = struct.Acc(:,3);
                fp = fft(PosTemp);
                fv = fft(VelTemp);
                Rpos = imag(fp);
                Rvel = imag(fv);

                Data_AC = [Data_AC; Rpos, Rvel, ones(length(Rpos),1)*gg];
            end
        end
    end
end

figure
gscatter(Data_AC(:,1),Data_AC(:,2),Data_AC(:,3))
xlabel('posac')
ylabel('velac')

DataOn = Data_AC(:,[1,2]);
labelac = Data_AC(:,[3]);

[Difference,ClassData,ProbData] = SimpleRelativeRBFTrain(DataOn,labelac);

figure
gscatter(DataOn(:,1),DataOn(:,2),labelac);
hold on
cc = 1;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),Difference{cc});
hold off
title('diffs from fft nov')


figure
gscatter(DataOn(:,1),DataOn(:,2),labelac);
hold on
cc = 2;
Surface3D(ClassData{cc}(:,1),ClassData{cc}(:,2),Difference{cc});
hold off
title('diffs from fft exp')


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
