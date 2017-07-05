% test classification using DPP on linear simulated data

% load up linear data

load LinearData11.mat

fsize = 14;

%% grab all data for plots

ratio = 0.3; %use half the data
coluse = [1,2,4];
mapping = [-1,1];

Data = [];
Labels = [];


for cc = 1:length(SegmentData.Class)
    for ii = 1:length(SegmentData.Class{cc}.Iteration)
        [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});
        
        %its too damn big
        ssmp = 1:100:ns;
        
        Data = [Data; SegmentData.Class{cc}.Iteration{ii}(ssmp,coluse)];
        Labels = [Labels; ones(length(ssmp),1)*mapping(cc)];
    end
end

figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
title('true class')

%% train time
tic
for t = 1:10
    t
    [~,subdata,sublabels] = GPRSubsample(Data,Labels, ratio);
    GPRMDL = GPRClassifyTrain(subdata,sublabels);
end
toc

%% online time
DataOn = Data(1:1000,:);
tic
for t = 1:10
    t
    [Class,ClassTime,ScoreTime]= GPRClassifyOnline(GPRMDL,DataOn);
end
toc
    
    
%% Get single class rbfs and plot

[Difference,ClassData,Model] = SimpleRelativeRBFTrain(Data,Labels);

%%Combine it
AllData = [ClassData{1}; ClassData{2}];
AllDiff = [Difference{1}; Difference{2}];

figure
gscatter(Data(:,1),Data(:,2),Labels)
hold on
Surface3D(AllData(:,1),AllData(:,2),AllDiff,'mesh');
hold off
xlabel('x','FontSize',fsize)
ylabel('U','FontSize',fsize)
zlabel('W_{RBF}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';
% title('DPP KL Divergence')
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'W_{RBF}','FontSize',fsize)

%% subsample that data

[subdiff,subdata,sublabels] = GPRSubsample(Data,Labels, ratio);

figure
p=gscatter(subdata(:,1),subdata(:,2),sublabels);

xlabel('x','FontSize',fsize)
ylabel('U','FontSize',fsize)
% zlabel('W_{KL}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';
xlim([0,4])
ylim([0,20])

%% train it

%Train model
GPRMDL = GPRClassifyTrain(subdata,sublabels);


%% sample classify

stashlabs = [];
stashscores = [];
stashdata = [];
for cc = 1:length(SegmentData.Class)
    for ii = 1:length(SegmentData.Class{cc}.Iteration)
        [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});
        ssmp = 1:100:ns;
        
        dtemp = SegmentData.Class{cc}.Iteration{ii}(ssmp,coluse);
        
        %now try a dummy classify
        [Class,ClassTime,ScoreTime]= GPRClassifyOnline(GPRMDL,dtemp);
        
        stashlabs = [ stashlabs; ClassTime];
        stashscores = [stashscores; ScoreTime];
        stashdata = [stashdata; dtemp];
    end
end


figure
gscatter(stashdata(:,1),stashdata(:,2),stashlabs,'rgc')
title('est class')


figure
scatter(stashdata(:,1),stashdata(:,2),10,stashscores)
xlabel('x','FontSize',fsize)
ylabel('U','FontSize',fsize)
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'L_{on}','FontSize',fsize)

%% perform leave one out

fsize = 14;

% get data
classes = length(SegmentData.Class);
runs = length(SegmentData.Class{1}.Iteration );
ratio = 0.3; %use half the data

%state and input columns
statecol = [1,2,4];
mapping = [-1,1];

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
ScoreTime = [];
DataAll = [];
ScoreAll = [];

%leave one out epoch
for oo = 1:runs
    
    fprintf('Epoch %d of %d \n',oo,runs)
    
    %clear our matrices
    Xtrain = [];
    Labelstrain = [];
    Xtest = [];
    Labelstest = [];
    
    %build up train and test data for this epoch
    for cc = 1:classes
        for ii = 1:runs
            [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});
            ssmp = 1:100:ns;
        
            tempstate = SegmentData.Class{cc}.Iteration{ii}(ssmp,statecol);
            
            
            if(ii == oo)
                %left out data
                Xtest = [Xtest; tempstate];
                Labelstest = [Labelstest; ones(length(ssmp),1)*mapping(cc) ];
            else
                %training data
                Xtrain = [Xtrain; tempstate];
                Labelstrain = [Labelstrain; ones(length(ssmp),1)*mapping(cc) ];
            end
        end
    end
    
    %subsample the data
    [subdiff,subdata,sublabels] = GPRSubsample(Xtrain,Labelstrain, ratio);
    
    %now we train for this epoch
    modelloo = GPRClassifyTrain(subdata,sublabels);

    %now we classify with test data
    for cc = 1:classes
        %grab one classes data
       xtemp = Xtest(Labelstest == mapping(cc),:);
       labeltrue = ones(size(xtemp,1),1)*mapping(cc);

       %try classifying 
       [cesttemp,cestcumtemp,scoretimetemp,salltemp] = GPRClassifyOnline(modelloo,xtemp);
        
       %find where the last time it changed classification was
       convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
       if(isempty(convtime))
           convtime = 0;
       end
       
       %store it temporarily
       ClassificationAll = [ClassificationAll; mapping(cc), cesttemp];
       ClassificationTimeAll = [ClassificationTimeAll; labeltrue, cestcumtemp];
       ConvergenceTime = [ConvergenceTime; convtime];
       DataAll = [DataAll; xtemp];
       ScoreTime.Class{cc}.run{oo} = scoretimetemp;
       ScoreAll = [ScoreAll; salltemp];
    end
        
end

avgconvtime = mean(ConvergenceTime)

corrtime = ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2);
acctime = mean(corrtime)
corr = ClassificationAll(:,1) == ClassificationAll(:,2);
acc = mean(corr)

%kappa=1.1 gets 72% time, 100%overall, 28 % convergence
%time
%kappa=1.2 gets 75% time, 100% classification, 36% convergence

%% Scores for correct and incorrect classifications

incorrectScores = abs(ScoreAll(ClassificationTimeAll(:,1) ~= ClassificationTimeAll(:,2) ));
correctScores = abs(ScoreAll(ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2) ));
C = [ones(length(incorrectScores),1)*1; ones(length(correctScores),1)*2] ;

mean(incorrectScores)
std(incorrectScores)
mean(correctScores)
std(correctScores)

figure
boxplot([incorrectScores; correctScores],C,'notch','on','labels',{'Incorrect','Correct'})
xlabel('Classification','FontSize',fsize)
ylabel('Weight \mu{*}','FontSize',fsize)


%% plot confidence over time too

figure
for cc = 1:classes
    for oo = 1:runs
        stemp = ScoreTime.Class{cc}.run{oo};
        nns = length(stemp);
        xper = ([1:nns] ./ nns )*100;
        if(mapping(cc) == -1)
            h1=plot(xper,stemp,'r-');
            hold on
        else
            h2=plot(xper,stemp,'c-');
            hold on
        end
    end
end
hold off
xlabel('% Trajectory','FontSize',fsize)
ylabel('L_{on}','FontSize',fsize)
lh=legend([h1(1),h2(1)],'-1','1');
lh.FontSize = fsize;
lh.Location = 'NorthWest';

