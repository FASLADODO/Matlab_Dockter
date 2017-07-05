%Segment the JigSaw Data

%this comes from loaddata.m
load('JIGSAWDATA.mat');

%% First lets just break it into Task, Skill Level, Surgeon, Trial

%clear it
FullData = [];

% for plotting skill level vs grs score
scoreStash = [];

%loop through all of Raw Data
for ii = 1:length(RawData.Task) 
    %need to update surgeon index for each skill level
    SkillSurgeonIDX = [0,0,0];
    for jj = 1:length(RawData.Task{ii}.Surgeon) 
        %skill level
        skilllevel = RawData.Task{ii}.Surgeon{jj}.SelfSkillNumber;
        %figure out current index
        SkillSurgeonIDX(skilllevel) = SkillSurgeonIDX(skilllevel) + 1;
        sid = SkillSurgeonIDX(skilllevel);
        %now through all trials from this surgeon
        for kk = 1:length(RawData.Task{ii}.Surgeon{jj}.Trial) 
            %current Data
            DataOn = RawData.Task{ii}.Surgeon{jj}.Trial{kk}.Data;
            InfoOn = RawData.Task{ii}.Surgeon{jj}.Trial{kk}.Info;
            
            %store it in fulldata
            FullData.Task{ii}.Skill{skilllevel}.Surgeon{sid}.Trial{kk}.Data = DataOn;
            FullData.Task{ii}.Skill{skilllevel}.Surgeon{sid}.Trial{kk}.Info = InfoOn;
        
            %for plots
            scoreStash = [scoreStash; InfoOn.All];
        end
    end
end

figure
gscatter(scoreStash(:,key.i.SelfProclaimedSkill),scoreStash(:,key.i.GRSSkillLevel))
xlabel('Self Proclaimed Skill')
ylabel('GRS Score')
title('score distirbution')

%% Now Lets do some segmentation

SegData = [];

%thresholds
thetathresh = -0.1;
velthresh = 0.0001;
minsteps = 10;



%loop through all data
for ii = 1:length(FullData.Task) 
    for jj = 1:length(FullData.Task{ii}.Skill) 
        for kk = 1:length(FullData.Task{ii}.Skill{jj}.Surgeon) 
            for ll = 1:length(FullData.Task{ii}.Skill{jj}.Surgeon{kk}) 
                %get the current data
                DataOn = FullData.Task{ii}.Skill{jj}.Surgeon{kk}.Trial{ll}.Data;
                InfoOn = FullData.Task{ii}.Skill{jj}.Surgeon{kk}.Trial{ll}.Info;
                SegData.Task{ii}.Skill{jj}.Surgeon{kk}.Trial{ll}.Info = InfoOn;
                [NN,SS] = size(DataOn);
                
                %indexes
                flipid = 0;
                startidx = 1;
                endidx = 1;
                segid = 1;
                
                %Left Hand
                velocities = DataOn(:,[key.c.Slave.L.dX,key.c.Slave.L.dY,key.c.Slave.L.dZ]);
                magvel = NormRowWise(velocities);
                gripperangle = DataOn(:,[key.c.Slave.L.GripAngle]);
                %find logical indices where data meets our specs
                isslow = magvel < velthresh;
                isclosed = gripperangle < thetathresh;
                isopened = gripperangle > thetathresh;
                isslowclose = (isslow + isclosed) > 1.5;
                isslowopen = (isslow + isopened) > 1.5;
                ischangetoclose = abs(diff(isslowclose));
                ischangetoopen = abs(diff(isslowopen));
                ischangeslow = abs(diff(isslow));
                
                figure(7)
                plotcolsL2 = [key.c.Time, key.c.Slave.L.GripAngle];
                scatter(DataOn(:,plotcolsL2(1)),DataOn(:,plotcolsL2(2)),'r.');
                DataStank = DataOn(ischangeslow==1,:);
                hold on
                scatter(DataStank(:,plotcolsL2(1)),DataStank(:,plotcolsL2(2)),'k+');
                hold off
       
                
                %loop through all time steps
                while tt < NN-1
                    tt = endidx + 1;
                    currentsearch = ischangeslow([tt:end],:);
                    endidx = find(currentsearch,1,'first') + tt;
                    if(isempty(endidx))
                       break; 
                    end
                    startidx = find(ischangeslow([startidx:endidx-1],:),1,'last') + startidx;
                    if(endidx - startidx > minsteps)
                        segmenttemp = DataOn([startidx:endidx],:);
                        SegData.Task{ii}.Skill{jj}.Surgeon{kk}.Trial{ll}.Segment{segid}.Data = segmenttemp;
                        segid = segid + 1;
                    end
                end
                
            end
        end
    end
end


            
            
            
            
            
