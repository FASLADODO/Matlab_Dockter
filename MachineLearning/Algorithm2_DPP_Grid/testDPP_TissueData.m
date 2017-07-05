%Test DPP grid on smart tool data

%load it up
load SmartToolSegments.mat


fsize = 14;

%% run through all grasps and store

fsize = 14;

DataAll = [];
Labels = [];
mapping = [-1 , 1];

%get grasp statistics

for ii = 1:5 %length(SegData.Donor)
    for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
        for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
            for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                nn = size(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data , 1);
                DataAll = [DataAll; SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data ];
                Labels = [Labels; ones(nn,1)*mapping(jj) ];
            end
        end
    end
end


%these are the columns we will use
pc = [key.c.Stress, key.c.Strain, key.c.dStrain]

%test columns
tc = [key.c.Stress, key.c.Strain, key.c.dStrain]
% tc = [key.c.Strain, key.c.dStrain]

figure
gscatter3(DataAll(:,tc(1)), DataAll(:,tc(2)), DataAll(:,tc(3)),Labels,'rc')
xlabel(key.c.all(tc(1)),'FontSize',fsize)
ylabel(key.c.all(tc(2)),'FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';

    
%% Now make the grid

coluse = [key.c.Stress, key.c.Strain];%, key.c.dStrain];

Data = DataAll(:,coluse);

colormapnew = flipud(cool);

split = 10;

[NN,SS] = size(Data);

[Model] = TrainDPPGrid(Data,Labels,split);

%% basic classify


stashlabs = [];
for ii = 1:5 %length(SegData.Donor)
    for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
        for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
            for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                dtemp = SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data(:,coluse);
                [Class,ClassTime,ScoreTime] = OnlineDPPGrid(dtemp,Model);
        
                stashlabs = [ stashlabs; ClassTime];
            end
        end
    end
end


figure
gscatter(Data(:,1),Data(:,2),stashlabs,'rgc')
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
xlabel('\theta','FontSize',fsize)
ylabel('F','FontSize',fsize)
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
xlabel('\theta','FontSize',fsize)
ylabel('F','FontSize',fsize)
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
xlabel('\theta','FontSize',fsize)
ylabel('F','FontSize',fsize)
zlabel('S_{j,k}','FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthWest';
%title('DPP KL Divergence')
colormap(flipud(cool))
hc = colorbar;
ylabel(hc, 'S_{j,k}','FontSize',fsize)


%% leave one patient out (inter patient variablitity)

% get data
classes = 2;
runs = 5;
split = 10;

%state and input functions (xdot + x + x^2
coluse = [key.c.Stress, key.c.Strain, key.c.dStrain];
mapping = [-1,1];

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
ScoreTime = [];
DataStore = [];
ScoreAll = [];

for oo = 1:runs %length(SegData.Donor)
    
    fprintf('Epoch %d of %d \n',oo,runs)
    
    %clear our matrices
    Xtrain = [];
    Labelstrain = [];
    Xtest = [];
    Labelstest = [];
    Ptrain = []; %parameters
    

    %stash just the training data
    for dd = 1:runs %length(SegData.Donor)
        for cc = 1:classes %tissue type
            for ll = 1:length(SegData.Donor{dd}.Tissue{cc}.Location) %locations
                for gg = length(SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp) %grasp
                    %get linearized version
                    tempstate = SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp{gg}.Data(:,coluse);

                    %number of points in this iteration
                    [ns,nc] = size(tempstate);


                    if(dd == oo)
                        %left out data
                        Xtest = [Xtest; tempstate];
                        Labelstest = [Labelstest; ones(ns,1)*mapping(cc) ];
                    else
                        %training data
                        Xtrain = [Xtrain; tempstate];
                        Labelstrain = [Labelstrain; ones(ns,1)*mapping(cc) ];
                    end
                end
            end
        end
    end
    
    %now we train for this epoch
    modelloo = TrainDPPGrid(Xtrain,Labelstrain,split);

    %now we classify with test data
    for cc = 1:classes

        %grab one classes data
       xtemp = Xtest(Labelstest ==  mapping(cc),:);
       labeltrue = ones(size(xtemp,1),1)*mapping(cc);

       %%try classifying 
       [cesttemp,cestcumtemp,scoretimetemp,salltemp] = OnlineDPPGrid(xtemp,modelloo);

       %find where the last time it changed classification was
       convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
       if(isempty(convtime) )
            convtime = 0;
       end
       %store it temporarily
       ClassificationAll = [ClassificationAll; mapping(cc), cesttemp];
       ClassificationTimeAll = [ClassificationTimeAll; labeltrue, cestcumtemp];
       ConvergenceTime = [ConvergenceTime; convtime, mapping(cc)];
       DataStore = [DataStore; xtemp];
       ScoreTime.Class{cc}.run{oo} = scoretimetemp;
       ScoreAll = [ScoreAll; salltemp];
    end
        
    
%     figure
%     gscatter(Xtrain(:,1),Xtrain(:,2),Labelstrain,'rc')
%     hold on
%     gscatter(Xtest(:,1),Xtest(:,2),Labelstest,'gk','ooo')
%     hold off
%     title('loo')

end

%convergence time
avgconvtime = mean(ConvergenceTime(:,1))
idxconv1 = find(ConvergenceTime(:,2) == -1);
idxconv2 = find(ConvergenceTime(:,2) == 1);
avgconvtime1 = mean(ConvergenceTime(idxconv1,1))
avgconvtime2 = mean(ConvergenceTime(idxconv2,1))

%inidividual accs (1: liver, 2 : pancreas)
idxTime1 = find(ClassificationTimeAll(:,1) == -1);
idxTime2 = find(ClassificationTimeAll(:,1) == 1);
corrtime1 = ClassificationTimeAll(idxTime1,1) == ClassificationTimeAll(idxTime1,2);
corrtime2 = ClassificationTimeAll(idxTime2,1) == ClassificationTimeAll(idxTime2,2);
acctime1 = mean(corrtime1)
acctime2 = mean(corrtime2)
idx1 = find(ClassificationAll(:,1) == -1);
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

%LODO gets 58 % time, 80% overall, 31% convergence

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

%% leave one Location out (intra-patient variability)

% get data
classes = 2;
runs = 5;
split = 10;

%state and input functions (xdot + x + x^2
coluse = [key.c.Stress, key.c.Strain, key.c.dStrain];
mapping = [-1,1];

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
ScoreTime = [];
DataStore = [];
ScoreAll = [];


for dd = 1:runs %length(SegData.Donor)
    
    fprintf('Epoch %d of %d \n',dd,runs)
    
    %now we leave out one location per
    sublocation = max([length(SegData.Donor{dd}.Tissue{1}.Location),length(SegData.Donor{dd}.Tissue{2}.Location)]);
    for oo = 1:sublocation
    
        %clear our matrices
        Xtrain = [];
        Labelstrain = [];
        Xtest = [];
        Labelstest = [];

        %stash just the training data
        for cc = 1:classes %tissue type
            %figure out which location we have to leave out
            %(different tissues have different # locations)
            numlocs = length(SegData.Donor{dd}.Tissue{cc}.Location);
            leaveout = mod(oo-1,numlocs)+1; %mod 1 indexed
            for ll = 1:numlocs %locations
                for gg = length(SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp) %grasp
                    %get linearized version
                    tempstate = SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp{gg}.Data(:,coluse);

                    %number of points in this iteration
                    [ns,nc] = size(tempstate);
                    
                    if(ll == leaveout)
                        %left out data
                        Xtest = [Xtest; tempstate];
                        Labelstest = [Labelstest; ones(ns,1)*mapping(cc) ];
                    else
                        %training data
                        Xtrain = [Xtrain; tempstate];
                        Labelstrain = [Labelstrain; ones(ns,1)*mapping(cc) ];
                    end
                end
            end
        end

        %now we train for this epoch
        modelloo = TrainDPPGrid(Xtrain,Labelstrain,split);

        %now we classify with test data
        for cc = 1:classes
           %grab one classes data
           xtemp = Xtest(Labelstest ==  mapping(cc),:);
           labeltrue = ones(size(xtemp,1),1)*mapping(cc);

           %try classifying
           [cesttemp,cestcumtemp,scoretimetemp,salltemp] = OnlineDPPGrid(xtemp,modelloo);

           %find where the last time it changed classification was
           convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
           if(isempty(convtime) )
               convtime = 0;
           end

           %store it temporarily
           ClassificationAll = [ClassificationAll; mapping(cc), cesttemp];
           ClassificationTimeAll = [ClassificationTimeAll; labeltrue, cestcumtemp];
           ConvergenceTime = [ConvergenceTime; convtime, mapping(cc)];
           DataStore = [DataStore; xtemp];
           ScoreTime.Class{cc}.run{oo} = scoretimetemp;
           ScoreAll = [ScoreAll; salltemp];
        end

    end
end

%convergence time
avgconvtime = mean(ConvergenceTime(:,1))
idxconv1 = find(ConvergenceTime(:,2) == -1);
idxconv2 = find(ConvergenceTime(:,2) == 1);
avgconvtime1 = mean(ConvergenceTime(idxconv1,1))
avgconvtime2 = mean(ConvergenceTime(idxconv2,1))

%inidividual accs (1: liver, 2 : pancreas)
idxTime1 = find(ClassificationTimeAll(:,1) == -1);
idxTime2 = find(ClassificationTimeAll(:,1) == 1);
corrtime1 = ClassificationTimeAll(idxTime1,1) == ClassificationTimeAll(idxTime1,2);
corrtime2 = ClassificationTimeAll(idxTime2,1) == ClassificationTimeAll(idxTime2,2);
acctime1 = mean(corrtime1)
acctime2 = mean(corrtime2)
idx1 = find(ClassificationAll(:,1) == -1);
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

%LOLO gets 59% time, 70% overall, 24% convergence time

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
        if(mapping(cc) == -1)
            h1=plot(stemp,'r-');
            hold on
        else
            h2=plot(stemp,'c-');
            hold on
        end
    end
end
hold off
xlabel('time (ms)','FontSize',fsize)
ylabel('M_{on}','FontSize',fsize)
lh=legend([h1(1),h2(1)],'-1','1');
lh.FontSize = fsize;
lh.Location = 'NorthWest';


