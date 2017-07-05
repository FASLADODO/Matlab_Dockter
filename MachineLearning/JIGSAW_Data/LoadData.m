clear key

%task type
key.t.knottying         = 1;
key.t.suturing          = 2;
key.t.needlepassing     = 3;
key.t.all   = {'Knot Tying','Suturing','Needle Passing'};

%skill level
key.s.Novice = 1;
key.s.Intermediate = 2;
key.s.Expert = 3;
key.s.all = {'Novice', 'Intermediate', 'Expert'};

%surgeon information
key.i.SelfProclaimedSkill   = 1;
key.i.GRSSkillLevel         = 2;
key.i.GRSRespectTissue      = 3;
key.i.GRSHandling           = 4;
key.i.GRSTimeMotion         = 5;
key.i.GRSFlow               = 6;
key.i.GRSOverall            = 7;
key.i.GRSQuality            = 8;
key.i.all = {'skill-level self-proclaimed', 'skill-level GRS', 'Respect for tissue',...
            'Suture/needle handling', 'Time and motion', 'Flow of operation',...
            'Overall performance', 'Quality of final product'};

%%% These are the raw values that come directly from the csv file %%%
%Master Left
key.c.Time      = 1; %This is extra
key.c.Master.L.X        = 2;
key.c.Master.L.Y        = 3;
key.c.Master.L.Z        = 4;
key.c.Master.L.Rot1     = 5;
key.c.Master.L.Rot2     = 6;
key.c.Master.L.Rot3     = 7;
key.c.Master.L.Rot4     = 8;
key.c.Master.L.Rot5     = 9;
key.c.Master.L.Rot6     = 10;
key.c.Master.L.Rot7     = 11;
key.c.Master.L.Rot8     = 12;
key.c.Master.L.Rot9     = 13;
key.c.Master.L.dX       = 14;
key.c.Master.L.dY       = 15;
key.c.Master.L.dZ       = 16;
key.c.Master.L.dRoll    = 17;
key.c.Master.L.dPitch   = 18;
key.c.Master.L.dYaw     = 19;
key.c.Master.L.GripAngle = 20;

%Master Right
key.c.Master.R.X        = 21;
key.c.Master.R.Y        = 22;
key.c.Master.R.Z        = 23;
key.c.Master.R.Rot1     = 24;
key.c.Master.R.Rot2     = 25;
key.c.Master.R.Rot3     = 26;
key.c.Master.R.Rot4     = 27;
key.c.Master.R.Rot5     = 28;
key.c.Master.R.Rot6     = 29;
key.c.Master.R.Rot7     = 30;
key.c.Master.R.Rot8     = 31;
key.c.Master.R.Rot9     = 32;
key.c.Master.R.dX       = 33;
key.c.Master.R.dY       = 34;
key.c.Master.R.dZ       = 35;
key.c.Master.R.dRoll    = 36;
key.c.Master.R.dPitch   = 37;
key.c.Master.R.dYaw     = 38;
key.c.Master.R.GripAngle = 39;

%Slave Left
key.c.Slave.L.X        = 40;
key.c.Slave.L.Y        = 41;
key.c.Slave.L.Z        = 42;
key.c.Slave.L.Rot1     = 43;
key.c.Slave.L.Rot2     = 44;
key.c.Slave.L.Rot3     = 45;
key.c.Slave.L.Rot4     = 46;
key.c.Slave.L.Rot5     = 47;
key.c.Slave.L.Rot6     = 48;
key.c.Slave.L.Rot7     = 49;
key.c.Slave.L.Rot8     = 50;
key.c.Slave.L.Rot9     = 51;
key.c.Slave.L.dX       = 52;
key.c.Slave.L.dY       = 53;
key.c.Slave.L.dZ       = 54;
key.c.Slave.L.dRoll    = 55;
key.c.Slave.L.dPitch   = 56;
key.c.Slave.L.dYaw     = 57;
key.c.Slave.L.GripAngle = 58;

%Slave Right
key.c.Slave.R.X        = 59;
key.c.Slave.R.Y        = 60;
key.c.Slave.R.Z        = 61;
key.c.Slave.R.Rot1     = 62;
key.c.Slave.R.Rot2     = 63;
key.c.Slave.R.Rot3     = 64;
key.c.Slave.R.Rot4     = 65;
key.c.Slave.R.Rot5     = 66;
key.c.Slave.R.Rot6     = 67;
key.c.Slave.R.Rot7     = 68;
key.c.Slave.R.Rot8     = 69;
key.c.Slave.R.Rot9     = 70;
key.c.Slave.R.dX       = 71;
key.c.Slave.R.dY       = 72;
key.c.Slave.R.dZ       = 73;
key.c.Slave.R.dRoll    = 74;
key.c.Slave.R.dPitch   = 75;
key.c.Slave.R.dYaw     = 76;
key.c.Slave.R.GripAngle = 77;


key.c.all = {'Time (us)', 'Master L X (mm)', 'Master L Y (mm)', 'Master L Z (mm)', 'Master L Rot1', 'Master L Rot2', 'Master L Rot3', 'Master L Rot4', 'Master L Rot5', 'Master L Rot6', 'Master L Rot7', 'Master L Rot8', 'Master L Rot9',...
           'Master L dX (mm/s)', 'Master L dY (mm/s)', 'Master L dZ (mm/s)', 'Master L dRoll', 'Master L dPitch', 'Master L dYaw', 'Master L Gripper Theta',...
           'Master R X (mm)', 'Master R Y (mm)', 'Master R Z (mm)', 'Master R Rot1', 'Master R Rot2', 'Master R Rot3', 'Master R Rot4', 'Master R Rot5', 'Master R Rot6', 'Master R Rot7', 'Master R Rot8', 'Master R Rot9',...
           'Master R dX (mm/s)', 'Master R dY (mm/s)', 'Master R dZ (mm/s)', 'Master R dRoll', 'Master R dPitch', 'Master R dYaw', 'Master R Gripper Theta',...
           'Slave L X (mm)', 'Slave L Y (mm)', 'Slave L Z (mm)', 'Slave L Rot1', 'Slave L Rot2', 'Slave L Rot3', 'Slave L Rot4', 'Slave L Rot5', 'Slave L Rot6', 'Slave L Rot7', 'Slave L Rot8', 'Slave L Rot9',...
           'Slave L dX (mm/s)', 'Slave L dY (mm/s)', 'Slave L dZ (mm/s)', 'Slave L dRoll', 'Slave L dPitch', 'Slave L dYaw', 'Slave L Gripper Theta',...
           'Slave R X (mm)', 'Slave R Y (mm)', 'Slave R Z (mm)', 'Slave R Rot1', 'Slave R Rot2', 'Slave R Rot3', 'Slave R Rot4', 'Slave R Rot5', 'Slave R Rot6', 'Slave R Rot7', 'Slave R Rot8', 'Slave R Rot9',...
           'Slave R dX (mm/s)', 'Slave R dY (mm/s)', 'Slave R dZ (mm/s)', 'Slave R dRoll', 'Slave R dPitch', 'Slave R dYaw', 'Slave R Gripper Theta'};

%% Load the data files from the csv files 

%first we direct it towards the three task folders within JIGSAW
directories = {'D:\temp\JIGSAW\Knot_Tying',...
                'D:\temp\JIGSAW\Suturing',...
                'D:\temp\JIGSAW\Needle_Passing'};
dataFolder = 'kinematics\AllGestures';
metafilestr = 'meta_file';
surgeonexpression = '[_](\w)\d{3}'; %for finding surgeon id in filename
levelstr = 'NIE';

RawData = [];
RawData.ReadMe = 'Data is stored by Task, Surgeon, Trial';
filecount = 0;

%loop through all the data directories
for jj=1:length(directories)
    %get all folders in this dir
    filesAndFolders = dir([directories{jj}]); 
    numOfFiles = length(filesAndFolders);
    
    %find the meta file
    metafileidx = 0;
    found = 0;
    for ff = 1:numOfFiles
          filename = filesAndFolders(ff).name;                               
          found = strfind(filename,metafilestr);         
          if ~isempty(found)
              fprintf('metafile: %s \n',filename);
              metafileidx = ff;
              break;  
          end                                       
    end
    
    %if we found a meta file now we access
    if(found)
        %load up the meta file
        tempname = [directories{jj} '\' filesAndFolders(metafileidx).name];
        fileID = fopen(tempname);
        metadata = textscan(fileID, '%s%s%f%f%f%f%f%f%f','delimiter','\t','MultipleDelimsAsOne',1);     
        fclose(fileID);

        %get the column with filenames
        datafilenames = metadata{1};
        numDataFiles = length(datafilenames);
        RawData.Task{jj}.MetaData = metadata;

        %surgeon index
        surgeonlist = '';
        surgeonidx = 0;
        trialidx = 1;
        
        %loop through filenames to get header info and raw data
        for ii = 1:numDataFiles
            %decide which file were doin now
            currentfile = datafilenames{ii};
            datafile = [directories{jj} '\' dataFolder '\' currentfile '.txt'];
            str = sprintf('reading file: %s.txt', currentfile);
            disp(str);
            
            %load up the data file
            DataTemp=load(datafile);
            %add the time data (30 hz)
            ln = size(DataTemp,1);
            timecol = [0:ln-1]';
            timecol = timecol.*(1/30);
            DataTemp = [timecol, DataTemp];
            
            %figure out which surgeon this is (identified by filename)
            [match] = regexp(currentfile,surgeonexpression,'tokens');
            surgeonalphabet = char(match{1,1});
            %increment surgeon ID if need be
            if(isempty(strfind(surgeonlist,surgeonalphabet)) )
                surgeonlist = strcat(surgeonlist,surgeonalphabet);
                surgeonidx = surgeonidx + 1;
                trialidx = 1;
            end
            
            %stash all of the data
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Data = DataTemp;
            
            %stash the surgeon info
            alphaskilllevel = char(metadata{2}(ii));
            numberskilllevel = strfind(levelstr,alphaskilllevel);
            infoall = [numberskilllevel,metadata{3}(ii),metadata{4}(ii),metadata{5}(ii),metadata{6}(ii),metadata{7}(ii),metadata{8}(ii),metadata{9}(ii)];
            %into the struct you go
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.FileName = currentfile;
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.All = infoall;
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.SkillLevelAlpha = alphaskilllevel;
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.SkillLevelNumber = numberskilllevel;
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.GRSLevel = metadata{3}(ii);
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.respect = metadata{4}(ii);
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.handling = metadata{5}(ii);
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.time = metadata{6}(ii);
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.flow = metadata{7}(ii);
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.overall = metadata{8}(ii);
            RawData.Task{jj}.Surgeon{surgeonidx}.Trial{trialidx}.Info.quality = metadata{9}(ii);
            RawData.Task{jj}.Surgeon{surgeonidx}.SelfSkillAlpha = alphaskilllevel;
            RawData.Task{jj}.Surgeon{surgeonidx}.SelfSkillNumber = numberskilllevel;
            RawData.Task{jj}.Surgeon{surgeonidx}.AlphaID = surgeonalphabet;
            
            %increment trial id
            trialidx = trialidx + 1;
            
            filecount = filecount + 1;
        end
    end
end

str = sprintf('Total Files = %d',filecount);
disp(str);

%% Test Plot

%get a random sample
ts = randi(length(RawData.Task),1);
ss = randi(length(RawData.Task{ts}.Surgeon),1);
rs = randi(length(RawData.Task{ts}.Surgeon{ss}.Trial),1);

%grab that data
DataOn = RawData.Task{ts}.Surgeon{ss}.Trial{rs}.Data;

%choose what to plot
plotcolsL = [key.c.Slave.L.X, key.c.Slave.L.Y, key.c.Slave.L.Z];
plotcolsR = [key.c.Slave.R.X, key.c.Slave.R.Y, key.c.Slave.R.Z];
plotcolsL2 = [key.c.Time, key.c.Slave.L.GripAngle];
plotcolsR2 = [key.c.Time, key.c.Slave.R.GripAngle];

%do it
figure
h1 = scatter3(DataOn(:,plotcolsL(1)),DataOn(:,plotcolsL(2)),DataOn(:,plotcolsL(3)),'r.');
hold on
h2 = scatter3(DataOn(:,plotcolsR(1)),DataOn(:,plotcolsR(2)),DataOn(:,plotcolsR(3)),'c.');
hold off
xlabel(key.c.all{plotcolsL(1)})
ylabel(key.c.all{plotcolsL(2)})
zlabel(key.c.all{plotcolsL(3)})
str = sprintf('Sample Plot, File: %s',RawData.Task{ts}.Surgeon{ss}.Trial{rs}.Info.FileName);
title(str)
legend([h1(1),h2(1)],'Left Hand','Right Hand')

figure
h1 = scatter(DataOn(:,plotcolsL2(1)),DataOn(:,plotcolsL2(2)),'r.');
hold on
h2 = scatter(DataOn(:,plotcolsR2(1)),DataOn(:,plotcolsR2(2)),'c.');
hold off
xlabel(key.c.all{plotcolsL2(1)})
ylabel(key.c.all{plotcolsL2(2)})
title('gripper angle')
legend([h1(1),h2(1)],'Left Hand','Right Hand')


%% Save the formatted data to a MATLAB binary file, for speed later.

save('JIGSAWDATA.mat','RawData','key')

