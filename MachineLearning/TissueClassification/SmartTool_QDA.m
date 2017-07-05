%% load up da data

load SmartToolSegments.mat
SegDataRaw = SegData;
%StateColumns = [key.c.dStrain,key.c.Strain,key.c.Stress];
StateColumns = [key.c.PosL,key.c.dPosL,key.c.ForceL];
plotcols = [key.c.PosL,key.c.dPosL,key.c.ForceL];
cslist = {key.t.all{1}; key.t.all{2}};
TotalDonors = 5; %ignore the last donor ("other")
TotalTissues = 2; %ignore the last tissue ("Open")


%% get overall QDA params

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
    
%train QDA
MdlQuadLoo = fitcdiscr(DataAll,LabelsAll,'DiscrimType','quadratic');

%estimate labels
LabelsEst = predict(MdlQuadLoo,DataAll);

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

pc2 = [plotcols(1), plotcols(3)]; 
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
                classest = predict(MdlQuadLoo,DataMat(:,StateColumns));
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


%% QDA classification with leave patient out (inter patient)

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
    clear ParamQuad
    ParamQuad = fitcdiscr(DataTrain,LabelsTrain,'DiscrimType','quadratic');
    
    %test classification with validation set
    for gg = 1:length(DataValidate.Tissue)
        for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
            dtemp = DataValidate.Tissue{gg}.Grasp{ss}(:,StateColumns);
            %predict using QDA
            classest = predict(ParamQuad,dtemp);
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

%% QDA classification with leave location out (intra patient)

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
    
    DonorAccuracy{pp} = [];
        
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
        clear ParamQuad
        ParamQuad = fitcdiscr(DataTrain,LabelsTrain,'DiscrimType','quadratic');

        tempacc = [];
        %test classification with validation set
        for gg = 1:length(DataValidate.Tissue)
            for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
                dtemp = DataValidate.Tissue{gg}.Grasp{ss}(:,StateColumns);
                %predict using QDA
                classest = predict(ParamQuad,dtemp);
                %overall estimate
                totalest = mean(classest);
                graspclass = (totalest > 1.5) + 1;
                GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
                classStore = [ classStore; classest, ones(length(classest),1)*gg];
                dataStore = [dataStore; DataValidate.Tissue{gg}.Grasp{ss}];
                
                tempacc = [tempacc;  graspclass, gg];
            end
        end    
        
        DonorAccuracy{pp} = [DonorAccuracy{pp}; tempacc];
    end
    
    corr = DonorAccuracy{pp}(:,1) == DonorAccuracy{pp}(:,2);
    acc = mean(corr);
    DonorClassification(pp) = acc;
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