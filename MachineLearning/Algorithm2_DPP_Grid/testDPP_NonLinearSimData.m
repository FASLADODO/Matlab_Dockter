% test classification using DPP on linear simulated data

% load up linear data

load NonLinearData12.mat

fsize = 14;

%% grab all data for plots

coluse = [1,4];
mapping = [-1,1];

Data = [];
Labels = [];


for cc = 1:length(SegmentData.Class)
    for ii = 1:length(SegmentData.Class{cc}.Iteration)
        [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});

        Data = [Data; SegmentData.Class{cc}.Iteration{ii}(:,coluse)];
        Labels = [Labels; ones(ns,1)*mapping(cc)];
    end
end

figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
title('true class')
    
%% Now make the grid

colormapnew = flipud(cool);

split = 13;

[NN,SS] = size(Data);

[Model] = TrainDPPGrid(Data,Labels,split);

%%
stashlabs = [];
for cc = 1:length(SegmentData.Class)
    for ii = 1:length(SegmentData.Class{cc}.Iteration)
        dtemp = SegmentData.Class{cc}.Iteration{ii}(:,coluse);
        %now try a dummy classify
        [Class,ClassTime,ScoreTime] = OnlineDPPGrid(dtemp,Model);
        
        stashlabs = [ stashlabs; ClassTime];
    end
end


figure
gscatter(Data(:,1),Data(:,2),stashlabs,'rc')
title('est class')
%% pretty plots


figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
% gscatter3(Data(:,1),Data(:,2),Data(:,3),Labels)
hold on
for ii = 1:size(Model.limits,1)
    minz = reshape(Model.limits(ii,1,:),1,SS);
    maxz = reshape(Model.limits(ii,2,:),1,SS);
    
    meanz = Model.means(ii,:);
    
    %get line order
    line([minz(1,1),minz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([maxz(1,1),maxz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [minz(1,2),minz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [maxz(1,2),maxz(1,2)])
    hold on
    %scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
xlabel('x','FontSize',fsize)
ylabel('U','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthWest';
% title('manual grid data')

%Try to plot KL
figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
limits = [min(Data); max(Data)]'
Surface3D(Model.means(:,1),Model.means(:,2),Model.KLRegion,'mesh',limits);
hold on
for ii = 1:size(Model.limits,1)
    minz = reshape(Model.limits(ii,1,:),1,SS);
    maxz = reshape(Model.limits(ii,2,:),1,SS);
    
    meanz = Model.means(ii,:);
    
    %get line order
    line([minz(1,1),minz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([maxz(1,1),maxz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [minz(1,2),minz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [maxz(1,2),maxz(1,2)])
    hold on
    %scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
xlabel('x','FontSize',fsize)
ylabel('U','FontSize',fsize)
zlabel('W_{KL}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';
% title('DPP KL Divergence')
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'W_{KL}','FontSize',fsize)



%compute S_{j,k} for kicks
% S = (Model.KLRegion .* Model.DensityRegion) ./ (abs(Model.KLRegion) + abs(Model.DensityRegion)) ;
figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
limits = [min(Data); max(Data)]'
Surface3D(Model.means(:,1),Model.means(:,2),Model.SWeightRegion,'mesh',limits);
hold on
for ii = 1:size(Model.limits,1)
    minz = reshape(Model.limits(ii,1,:),1,SS);
    maxz = reshape(Model.limits(ii,2,:),1,SS);
    
    meanz = Model.means(ii,:);
    
    %get line order
    line([minz(1,1),minz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([maxz(1,1),maxz(1,1)], [minz(1,2),maxz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [minz(1,2),minz(1,2)])
    hold on
    line([minz(1,1),maxz(1,1)], [maxz(1,2),maxz(1,2)])
    hold on
    %scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
xlabel('x','FontSize',fsize)
ylabel('U','FontSize',fsize)
zlabel('S_{j,k}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthWest';
%title('DPP KL Divergence')
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'S_{j,k}','FontSize',fsize)



%% perform leave one out

fsize = 14;

% get data
classes = length(SegmentData.Class);
runs = length(SegmentData.Class{1}.Iteration );
split = 13;

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
    Ptrain = []; %parameters
    
    %build up train and test data for this epoch
    for cc = 1:classes
        for ii = 1:runs
            %number of points in this iteration
            [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});
            
            if(ii == oo)
                %left out data
                Xtest = [Xtest; SegmentData.Class{cc}.Iteration{ii}(:,statecol)];
                Labelstest = [Labelstest; ones(ns,1)*mapping(cc) ];
            else
                %training data
                Xtrain = [Xtrain; SegmentData.Class{cc}.Iteration{ii}(:,statecol)];
                Labelstrain = [Labelstrain; ones(ns,1)*mapping(cc) ];
            end
        end
    end
    


    %now we train for this epoch
    modelloo = TrainDPPGrid(Xtrain,Labelstrain,split);
    
    %now we classify with test data
    for cc = 1:classes
        %grab one classes data
       xtemp = Xtest(Labelstest == mapping(cc),:);
       labeltrue = ones(size(xtemp,1),1)*mapping(cc);

       %try classifying 
       [cesttemp,cestcumtemp,scoretimetemp,salltemp] = OnlineDPPGrid(xtemp,modelloo);
        
       %find where the last time it changed classification was
       convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);

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

%kappa=1.1 gets 78.6% time, 100%overall, 25 % convergence
%time
%kappa=1.2 gets 82.8 % time, 100% overall, 19% convergence


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
ylabel('Weight S_{j,k}','FontSize',fsize)


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
ylabel('M_{on}','FontSize',fsize)
lh=legend([h1(1),h2(1)],'-1','1');
lh.FontSize = fsize;
lh.Location = 'NorthWest';

