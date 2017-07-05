%% load up da data

load SmartToolSegments.mat
SegDataRaw = SegData;
% StateColumns = [key.c.dStrain,key.c.StrainFilt,key.c.ForceLFilt,key.c.ResistFilt];
StateColumns = [key.c.PosL,key.c.dPosL,key.c.ForceL];
plotcols = [key.c.PosL,key.c.dPosL,key.c.ForceL];
cslist = {key.t.all{1}; key.t.all{2}};
TotalDonors = 5; %ignore the last donor ("other")
TotalTissues = 2; %ignore the last tissue ("Open")
rounds = 20; %boosting rounds



%% get overall random forest params

DataAll = [];
LabelsAll = [];
DataPlot =[];

for pp = 1:TotalDonors
    for gg = 1:TotalTissues
        for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
            for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                DataMat = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Data;

                DataAll = [DataAll; DataMat(:,StateColumns)];
                LabelsAll = [LabelsAll; ones(size(DataMat,1),1)*gg];
                
                DataPlot = [DataPlot; DataMat];
            end
        end
    end
end
    
%train random forest
tic
MdlAll = TreeBagger(rounds,DataAll,LabelsAll)
toc


%estimate labels
%estimate labels using the forest
tic
Est = predict(MdlAll,DataAll); 
LabelsEst = cellfun(@str2double, Est);
toc

figure
gscatter3(DataPlot(:,plotcols(1)),DataPlot(:,plotcols(2)),DataPlot(:,plotcols(3)),LabelsAll)
xlabel(key.c.all(plotcols(1)))
ylabel(key.c.all(plotcols(2)))
ylabel(key.c.all(plotcols(3)))
title('true class')

figure
gscatter3(DataPlot(:,plotcols(1)),DataPlot(:,plotcols(2)),DataPlot(:,plotcols(3)),LabelsEst)
xlabel(key.c.all(plotcols(1)))
ylabel(key.c.all(plotcols(2)))
ylabel(key.c.all(plotcols(3)))
title('Est class')

%% plot each donor seperately

pc2 = [key.c.PosL,key.c.ForceL];
% pc2 = [key.c.Strain,key.c.Stress];

% get all data from all grasps as a test
for pp = 1:TotalDonors
    DataDonor{pp} = [];
    LabelDonor{pp} = [];
    LabelEstDonor{pp} = [];
    for gg = 1:TotalTissues
        for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
            for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                DataMat = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Data;

                DataDonor{pp} = [DataDonor{pp}; DataMat];
                LabelDonor{pp} = [LabelDonor{pp}; ones(size(DataMat,1),1)*gg];
                
                %estimate labels
                classest = predict(MdlAll,DataMat(:,StateColumns));
                LabelEstDonor{pp} = [LabelEstDonor{pp}; classest];
            end
        end
    end
    

    figure
    gscatter(DataDonor{pp}(:,pc2(1)),DataDonor{pp}(:,pc2(2)),LabelDonor{pp})
    xlabel(key.c.all(pc2(1)))
    ylabel(key.c.all(pc2(2)))
    str = sprintf('%s, true class',key.d.all{pp});
    title(str)
    
    figure
    gscatter(DataDonor{pp}(:,pc2(1)),DataDonor{pp}(:,pc2(2)),LabelEstDonor{pp})
    xlabel(key.c.all(pc2(1)))
    ylabel(key.c.all(pc2(2)))
    str = sprintf('%s, est class',key.d.all{pp});
    title(str)
end


%% forest classification with leave patient out (inter patient)

epochs = TotalDonors;
trainIDX = 1:epochs;
valIDX = 1;

%storage
classStore = [];
GraspClassStore = [];
dataStore = [];

for oo = 1:epochs
    %decide on our leave one out
    trainIDX = 1:epochs;
    trainIDX(oo) = [];
    valIDX = oo;
    
    %for storage
    DataTrain = [];
    LabelsTrain = [];
    DataValidate = [];
    LabelsValidate = [];
    
    
    %get training and validation data
    for pp = 1:TotalDonors
        for gg = 1:TotalTissues
            segid = 1;
            for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
                for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                    DataMat = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Data;
                    
                    if(any(trainIDX==pp))
                        DataTrain = [DataTrain; DataMat(:,StateColumns)];
                        LabelsTrain = [LabelsTrain; ones(size(DataMat,1),1)*gg ];
                    elseif(any(valIDX==pp))
                        DataValidate.Tissue{gg}.Grasp{segid} = [DataMat];
                        LabelsValidate.Tissue{gg}.Grasp{segid} = [ones(size(DataMat,1),1)*gg ];
                        segid = segid + 1;
                    end
                end
            end
        end
    end

    %Train params with training set
    clear ForestModel
    ForestModel = TreeBagger(rounds,DataTrain,LabelsTrain);

    %test classification with validation set
    for gg = 1:length(DataValidate.Tissue)
        for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
            dtemp = DataValidate.Tissue{gg}.Grasp{ss}(:,StateColumns);
            %predict using random forest
            Est = predict(ForestModel,dtemp); 
            classest = cellfun(@str2double, Est);
            %overall estimate
            totalest = mean(classest);
            graspclass = (totalest > 1.5) + 1;
            GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
            classStore = [ classStore; classest, ones(length(classest),1)*gg];
            dataStore = [dataStore; DataValidate.Tissue{gg}.Grasp{ss}];
        end
    end
end


corrall = classStore(:,1) == classStore(:,2);
accall = mean(corrall)

corrgrasp = GraspClassStore(:,1) == GraspClassStore(:,2);
accgrasp = mean(corrgrasp)

figure
gscatter3(dataStore(:,plotcols(1)),dataStore(:,plotcols(2)),dataStore(:,plotcols(3)),classStore(:,1))
xlabel(key.c.all(plotcols(1)))
ylabel(key.c.all(plotcols(2)))
ylabel(key.c.all(plotcols(3)))
title('Est Class Inter')

%% forest classification with leave location out (intra patient)

epochs = length(SegData.Donor{1}.Tissue{1}.Location);
trainIDX = 1:epochs;
valIDX = 1;

%storage
errorStore = [];
classStore = [];
GraspClassStore = [];
dataStore = [];

%get training and validation data
for pp = 1:TotalDonors
    
    %decide on our leave one out based on largest number of locations
    for gg = 1:TotalTissues
        epochall(gg) = length(SegData.Donor{pp}.Tissue{gg}.Location );
    end
    epochs=max(epochall);
        
    for oo = 1:epochs
        %for storage
        DataTrain = [];
        LabelsTrain = [];
        DataValidate = [];
        LabelsValidate = [];
        
        %loop through all tissues, locations, and grasps
        for gg = 1:TotalTissues
            %the number of locations left out may change
            dropidx = mod(oo-1,epochall(gg))+1; %because some tissues may have fewer locations
            trainIDX = 1:epochall(gg);
            trainIDX(dropidx) = [];
            valIDX = dropidx;
            segid = 1;
            for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
                for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                    %grab the current data
                    DataMat = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Data;
                    
                    if(any(trainIDX==ll))
                        DataTrain = [DataTrain; DataMat(:,StateColumns)];
                        LabelsTrain = [LabelsTrain; ones(size(DataMat,1),1)*gg ];
                    elseif(any(valIDX==ll))
                        DataValidate.Tissue{gg}.Grasp{segid} = [DataMat];
                        LabelsValidate.Tissue{gg}.Grasp{segid} = [ones(size(DataMat,1),1)*gg ];
                        segid = segid + 1;
                    end
                end
            end
        end


        %Train params with training set
        clear ForestModel
        ForestModel = TreeBagger(rounds,DataTrain,LabelsTrain);

        %test classification with validation set
        for gg = 1:length(DataValidate.Tissue)
            for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
                dtemp = DataValidate.Tissue{gg}.Grasp{ss}(:,StateColumns);
                %predict using random forest
                Est = predict(ForestModel,dtemp); 
                classest = cellfun(@str2double, Est);
                %overall estimate
                totalest = mean(classest);
                graspclass = (totalest > 1.5) + 1;
                GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
                classStore = [ classStore; classest, ones(length(classest),1)*gg];
                dataStore = [dataStore; DataValidate.Tissue{gg}.Grasp{ss}];
            end
        end       
    end
end


corrall = classStore(:,1) == classStore(:,2);
accall = mean(corrall)

corrgrasp = GraspClassStore(:,1) == GraspClassStore(:,2);
accgrasp = mean(corrgrasp)

figure
gscatter3(dataStore(:,plotcols(1)),dataStore(:,plotcols(2)),dataStore(:,plotcols(3)),classStore(:,1))
xlabel(key.c.all(plotcols(1)))
ylabel(key.c.all(plotcols(2)))
ylabel(key.c.all(plotcols(3)))
title('Est Class Intra')