%Test DPP grid on smart tool data

%load it up
load SmartToolSegments.mat

%% get grasp statistics

TotalLocations = zeros(1,2);
TotalGrasps = zeros(1,2);
GraspLengths = [];

id1 = 1;
id2 = 1;

for ii = 1:5 %length(SegData.Donor)
    for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
        TotalLocations(jj) = TotalLocations(jj) + length(SegData.Donor{ii}.Tissue{jj}.Location);
        for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
            TotalGrasps(jj) = TotalGrasps(jj) + length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp);
            for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                temp = SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data;
                nn = size(temp , 1);
                if(jj == 1)
                    GraspLengths(id1,jj) = nn;
                    id1 = id1 + 1;
                elseif(jj == 2)
                    GraspLengths(id2,jj) = nn;
                    id2 = id2 + 1;
                end
                
            end
        end
    end
end
TotalLocations
TotalGrasps
mean(GraspLengths)
std(GraspLengths)
mean(GraspLengths(:))
std(GraspLengths(:))

%% run through all grasps and store

DataAll = [];
Labels = [];
map = [-1,1];

%get grasp statistics

for ii = 1:5 %length(SegData.Donor)
    for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
        for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
            for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                nn = size(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data , 1);
                DataAll = [DataAll; SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data ];
                Labels = [Labels; ones(nn,1)*map(jj) ];
            end
        end
    end
end


%these are the columns we will use
pc = [key.c.Stress, key.c.Strain]

%test columns
tc = [key.c.Stress, key.c.Strain]
% tc = [key.c.Strain, key.c.dStrain]
figure
gscatter(DataAll(:,tc(1)), DataAll(:,tc(2)),Labels,'rc')
xlabel(key.c.all(tc(1)))
ylabel(key.c.all(tc(2)))


%% Train gridded weak 

split = 11;


Data = DataAll(:,pc);
[NN,SS] = size(Data);

[Model] = TrainDPPGridWeak(Data,Labels,split);

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
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('manual grid data')

%Try to plot surface with density
%plot grid on data
figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
Surface3D(Model.means(:,1),Model.means(:,2),mean(Model.WeightRegion,2));
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
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('DPP Weights')

%% Try it online (weak)

[ClassEst,RawStore] = OnlineDPPGridWeak(Data,Model);

figure
gscatter(Data(:,1),Data(:,2),ClassEst(:,1),'rgc')
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
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('Class Est')

%% plot raw sums (weak)

figure
scatter(Data(:,1),Data(:,2),10,RawStore)
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
end
hold off
title('Class Weights')
colorbar 
colormap cool

%% leave one out (weak)

DataLoo = [];
ValLoo = [];

for oo = 1:5 %length(SegData.Donor)
    DataTrain = [];
    LabelsTrain = [];

    %stash just the training data
    for ii = 1:5 %length(SegData.Donor)
        if(ii ~= oo)
            for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
                for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
                    for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                        dtemp = SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data(:,pc);
                        nn = size(dtemp , 1);
                        
                        DataTrain = [DataTrain; dtemp];
                        LabelsTrain = [LabelsTrain; ones(nn,1)*map(jj)];
                    end
                end
            end
        end
    end
    
    %train the model
    [ModelLoo] = TrainDPPGridWeak(DataTrain,LabelsTrain,split);

    DataTempLoo = [];
    LabelsTempLoo = [];
    %do leave on out
    for ii = oo
        for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
            for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
                for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                    dtemp = SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data(:,pc);
                    nn = size(dtemp , 1);
                    
                    %classify with the online dpp
                    [SampleClass,SampleScore] = OnlineDPPGridWeak(dtemp,ModelLoo);
                    
                    %Overall Grap Estimate
                    GrapClass = Model.f_classify(mean(SampleScore));
                    
                    %cum sum weighted score
                    GraspScore = cumsum(SampleScore);
                    
                    %store 
                    DataLoo = [DataLoo; dtemp];
                    ValLoo = [ValLoo; ones(nn,1)*map(jj), SampleClass, ones(nn,1)*GrapClass, SampleScore, GraspScore ];
                    
                    %just for plots
                    DataTempLoo = [DataTempLoo; dtemp];
                    LabelsTempLoo = [LabelsTempLoo; ones(nn,1)*map(jj)];
                end
            end
        end
    end
    
    figure
    gscatter(DataTrain(:,1),DataTrain(:,2),LabelsTrain,'rc')
    hold on
    gscatter(DataTempLoo(:,1),DataTempLoo(:,2),LabelsTempLoo,'gk','ooo')
    hold off
    title('loo')

end

figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
title('True Class')

figure
gscatter(DataLoo(:,1),DataLoo(:,2),ValLoo(:,2),'rgc')
title('Est. class (per sample)')

figure
gscatter(DataLoo(:,1),DataLoo(:,2),ValLoo(:,3),'rc')
title('Est. class (per grasp)')

figure
scatter(DataLoo(:,1),DataLoo(:,2),10,ValLoo(:,4))
title('Weighted Score (per sample)')
colorbar
colormap cool

figure
scatter(DataLoo(:,1),DataLoo(:,2),10,ValLoo(:,5))
title('Weighted Score (per grasp)')
colorbar
colormap cool


samplecorr = ValLoo(:,1) == ValLoo(:,2);
sampleacc = mean(samplecorr)

graspcorr = ValLoo(:,1) == ValLoo(:,3);
graspacc = mean(graspcorr)

% Should give ~ 62 %
return;

%% Trained Gridded Prob

split = 11;

Data = DataAll(:,pc);
[NN,SS] = size(Data);

%weaponizes the training
[Model] = TrainDPPGrid(Data,Labels,split);

%plot grid on data
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
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('manual grid data')


%Try to plot surface with density
%plot grid on data
figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
hold on
Surface3D(Model.means(:,1),Model.means(:,2),Model.DensityRegion);
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
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('DPP Density')

%Try to plot surface with probs
%plot grid on data
figure
gscatter(Data(:,1),Data(:,2),Labels)
hold on
Surface3D(Model.means(:,1),Model.means(:,2),Model.KLRegion(:,1));
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
    scatter(meanz(1),meanz(2),'k+')
    hold on
end
hold off
title('DPP KL')

%% Try to classify each grasp online all at once (prob)

DataOnline = [];
ValOnline = [];

for ii = 1:5 %length(SegData.Donor)
    for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
        for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
            for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                dtemp = SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data(:,pc);
                nn = size(dtemp , 1);
               
                %do online stuff
                [Class,Score,ScoreTime] = OnlineDPPGrid(dtemp,Model);
                
                %store 
                DataOnline = [DataOnline; dtemp];
                ValOnline = [ValOnline; ones(nn,1)*Class, ones(nn,1)*jj, ones(nn,1)*Score];
            end
        end
    end
end

figure
gscatter(Data(:,1),Data(:,2),Labels)
title('true class')

figure
gscatter(DataOnline(:,1),DataOnline(:,2),ValOnline(:,1))
title('Est class')

figure
scatter(DataOnline(:,1),DataOnline(:,2),10,ValOnline(:,3))
title('final val')
colorbar
colormap cool

corr = ValOnline(:,1) == ValOnline(:,2);
acc = mean(corr)


%% leave one out

DataLoo = [];
ValLoo = [];

for oo = 1:5 %length(SegData.Donor)
    DataTrain = [];
    LabelsTrain = [];

    %stash just the training data
    for ii = 1:5 %length(SegData.Donor)
        if(ii ~= oo)
            for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
                for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
                    for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                        dtemp = SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data(:,pc);
                        nn = size(dtemp , 1);
                        
                        DataTrain = [DataTrain; dtemp];
                        LabelsTrain = [LabelsTrain; ones(nn,1)*jj];
                    end
                end
            end
        end
    end
    
    %train the model
    [ModelLoo] = TrainDPPGrid(DataTrain,LabelsTrain,split);

    DataTempLoo = [];
    LabelsTempLoo = [];
    %do leave on out
    for ii = oo
        for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
            for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
                for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                    dtemp = SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data(:,pc);
                    nn = size(dtemp , 1);
                    
                    [Class,Score,ScoreTime] = OnlineDPPGrid(dtemp,ModelLoo);

                    %store 
                    DataLoo = [DataLoo; dtemp];
                    ValLoo = [ValLoo; ones(nn,1)*Class, ones(nn,1)*jj, ones(nn,1)*Score, ScoreTime];
                    
                    %just for plots
                    DataTempLoo = [DataTempLoo; dtemp];
                    LabelsTempLoo = [LabelsTempLoo; ones(nn,1)*jj];
                end
            end
        end
    end
    
    figure
    gscatter(DataTrain(:,1),DataTrain(:,2),LabelsTrain)
    hold on
    gscatter(DataTempLoo(:,1),DataTempLoo(:,2),LabelsTempLoo,'kg','oo')
    hold off
    title('loo')

end

figure
gscatter(Data(:,1),Data(:,2),Labels)
title('true class')

figure
gscatter(DataLoo(:,1),DataLoo(:,2),ValLoo(:,1))
title('Est class')

figure
scatter(DataLoo(:,1),DataLoo(:,2),10,ValLoo(:,4))
title('Time Val')
colorbar
colormap cool

figure
scatter(DataLoo(:,1),DataLoo(:,2),10,ValLoo(:,3))
title('final val')
colorbar
colormap cool

corr = ValLoo(:,1) == ValLoo(:,2);
acc = mean(corr)

