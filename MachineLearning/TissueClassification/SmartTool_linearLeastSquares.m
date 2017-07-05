%% load up da data

load SmartToolSegments.mat
SegDataRaw = SegData;
StateColumns = [key.c.dPosL,key.c.PosL,key.c.PosL2,key.c.PosL4];
% StateColumns = [key.c.dStrain,key.c.Strain,key.c.Strain2,key.c.Strain4];
InputColumns = [key.c.Stress];
% plotcols = [key.c.Strain,key.c.dStrain,key.c.Stress];
plotcols = [key.c.PosL,key.c.dPosL,key.c.ForceL];
cslist = {key.t.all{1}; key.t.all{2}};
TotalDonors = 5; %ignore the last donor ("other")
TotalTissues = 2; %ignore the last tissue ("Open")


%% get true linear params


for gg = 1:TotalTissues
    ClassData{gg} = [];
    ClassInput{gg} = [];
    AllData{gg} = [];
end
% get all data from all grasps as a test
for pp = 1:TotalDonors
    for gg = 1:TotalTissues
        for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
            for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                DataMat = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Data;

                ClassData{gg} = [ClassData{gg}; DataMat(:,StateColumns)];
                ClassInput{gg} = [ ClassInput{gg}; DataMat(:,InputColumns)];
                
                AllData{gg} = [AllData{gg}; DataMat];
            end
        end
    end
end

%run with true params for all data
for gg = 1:TotalTissues
    dtemp = ClassData{gg};
    utemp = ClassInput{gg};
    phiTrue{gg} = pinv(dtemp)*utemp;
end

figure
gg = 1;
scatter3(AllData{gg}(:,plotcols(1)),AllData{gg}(:,plotcols(2)),AllData{gg}(:,plotcols(3)),'b.')
hold on
gg = 2;
scatter3(AllData{gg}(:,plotcols(1)),AllData{gg}(:,plotcols(2)),AllData{gg}(:,plotcols(3)),'r.')
hold off
xlabel(key.c.all(plotcols(1)))
ylabel(key.c.all(plotcols(2)))
zlabel(key.c.all(plotcols(3)))
title('all data true class')

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
                
                err = [];
                %classify it too
                for gc = 1:TotalTissues
                    dtemp = DataMat(:,StateColumns);
                    utemp = DataMat(:,InputColumns);
                    errtemp = abs(utemp - dtemp*phiTrue{gc});
                    err = [err, errtemp];
                end
                [~,classest] = min(err,[],2);
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


%% least squares classification with leave patient out (inter patient)

epochs = TotalDonors;
trainIDX = 1:epochs;
runIDX = 1;

%storage
errorStore = [];
classStore = [];
GraspClassStore = [];
dataStore = [];

for oo = 1:epochs
    %decide on our leave one out
    trainIDX = 1:epochs;
    trainIDX(oo) = [];
    valIDX = oo;
    
    %for storage
    for gg = 1:TotalTissues
        DataTrain{gg} = [];
        InputTrain{gg} = [];
        LabelsTrain{gg} = [];
    end
    DataValidate = [];
    InputValidate = [];
    LabelsValidate = [];
    
    
    %get training and validation data
    for pp = 1:TotalDonors
        for gg = 1:TotalTissues
            segid = 1;
            for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
                for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                    DataMat = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss}.Data;
                    
                    if(any(trainIDX==pp))
                        DataTrain{gg} = [DataTrain{gg}; DataMat(:,StateColumns)];
                        InputTrain{gg} = [ InputTrain{gg}; DataMat(:,InputColumns)];
                    elseif(any(valIDX==pp))
                        DataValidate.Tissue{gg}.Grasp{segid} = DataMat;
                        InputValidate.Tissue{gg}.Grasp{segid} = DataMat;
                        segid = segid + 1;
                    end
                end
            end
        end
    end

    %Train params with training set
    for gg = 1:TotalTissues
        dtemp =  DataTrain{gg};
        utemp = InputTrain{gg};
        ParamsTrain{gg} = pinv(dtemp)*utemp;
    end
    
    %test classification with validation set
    for gg = 1:length(DataValidate.Tissue)
        for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
            dtemp = DataValidate.Tissue{gg}.Grasp{ss}(:,StateColumns);
            utemp = InputValidate.Tissue{gg}.Grasp{ss}(:,InputColumns);
            error = [];
            for cc = 1:TotalTissues
                error = [error, abs(utemp - dtemp*ParamsTrain{cc})];
            end
            errorStore = [errorStore; error];
            [~,classest] = min(error,[],2);
            totalest = mean(classest);
            graspclass = (totalest > 1.5) + 1;
            GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
            classStore = [ classStore; classest, ones(length(error),1)*gg];
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
zlabel(key.c.all(plotcols(3)))
title('Est Class Inter')

%% least squares classification with leave location out (intra patient)

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
        for gg = 1:TotalTissues
            DataTrain{gg} = [];
            InputTrain{gg} = [];
            LabelsTrain{gg} = [];
        end
        DataValidate = [];
        InputValidate = [];
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
                        DataTrain{gg} = [DataTrain{gg}; DataMat(:,StateColumns)];
                        InputTrain{gg} = [ InputTrain{gg}; DataMat(:,InputColumns)];
                    elseif(any(valIDX==ll))
                        DataValidate.Tissue{gg}.Grasp{segid} = DataMat;
                        InputValidate.Tissue{gg}.Grasp{segid} = DataMat;
                        segid = segid + 1;
                    end
                end
            end
        end


        %Train params with training set
        for gg = 1:TotalTissues
            dtemp =  DataTrain{gg};
            utemp = InputTrain{gg};
            %least squares parameters
            ParamsTrain{gg} = pinv(dtemp)*utemp;
        end
        
        tempacc = [];

        %test classification with validation set
        for gg = 1:length(DataValidate.Tissue)
            for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
                dtemp = DataValidate.Tissue{gg}.Grasp{ss}(:,StateColumns);
                utemp = InputValidate.Tissue{gg}.Grasp{ss}(:,InputColumns);
                %classify using least squares
                error = [];
                for cc = 1:TotalTissues
                    error = [error, abs(utemp - dtemp*ParamsTrain{cc})];
                end
                errorStore = [errorStore; error];
                [~,classest] = min(error,[],2);
                totalest = mean(classest);
                graspclass = (totalest > 1.5) + 1;
                GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
                classStore = [ classStore; classest, ones(length(error),1)*gg];
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
zlabel(key.c.all(plotcols(3)))
title('Est Class Intra')

%% least squares classification intra location

clear DataTrain InputTrain DataValidate InputValidate

epochs = length(SegData.Donor{1}.Tissue{1}.Location{1}.Grasp);
trainIDX = 1:epochs;
valIDX = 1;

%storage
errorStore = [];
classStore = [];
GraspClassStore = [];
dataStore = [];

totallocs = 5;

allAll = [];

%get training and validation data
for pp = 1:TotalDonors

    stashLocations = [];
    stashLabels = [];
    accDonor = [];
    %loop through all tissues, locations, and grasps
    for ll = 1:totallocs

         for oo = 1:epochs

            %for storage
            for ggc = 1:TotalTissues
                DataTrain{ggc} = [];
                InputTrain{ggc} = [];
                DataValidate{ggc} = [];
                InputValidate{ggc} = [];

            end
            %what we leave out
            trainIDX = 1:epochs;
            trainIDX(oo) = [];
            valIDX = oo;


            for gg = 1:TotalTissues
                
                locid = ll; %because not all tissues have 5 locations
                if(locid > length(SegData.Donor{pp}.Tissue{gg}.Location) )
                    locid = mod(ll-1,length(SegData.Donor{pp}.Tissue{gg}.Location))+1;
                end

                for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{locid}.Grasp )
                    %we have to regrab all tissues


                    %grab the current data
                    DataMat = SegData.Donor{pp}.Tissue{gg}.Location{locid}.Grasp{ss}.Data;

                    if(any(trainIDX==ss))
                        DataTrain{gg} = [ DataMat(:,StateColumns)];
                        InputTrain{gg} = [ DataMat(:,InputColumns)];
                    elseif(any(valIDX==ss))
                        DataValidate{gg} = DataMat;
                        InputValidate{gg} = DataMat;
                    end
                end
            end


            %Train params with training set
            for gg = 1:TotalTissues
                dtemp =  DataTrain{gg};
                utemp = InputTrain{gg};
                %least squares parameters
                ParamsTrain{gg} = pinv(dtemp)*utemp;
            end

            for gg = 1:TotalTissues
                %test classification with validation set
                dtemp = DataValidate{gg}(:,StateColumns);
                utemp = InputValidate{gg}(:,InputColumns);
                %classify using least squares
                error = [];
                for cc = 1:TotalTissues
                    error = [error, abs(utemp - dtemp*ParamsTrain{cc})];
                end
                errorStore = [errorStore; error];
                [~,classest] = min(error,[],2);
                totalest = mean(classest);
                graspclass = (totalest > 1.5) + 1;
                GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
                classStore = [ classStore; classest, ones(length(error),1)*gg];
                dataStore = [dataStore; DataValidate{gg}];

                stashLocations = [stashLocations; DataValidate{gg}];
                stashLabels = [stashLabels; classest];

                accDonor = [accDonor;  graspclass, gg];
            end
         end
    end
    
    corrt = accDonor(:,1) == accDonor(:,2);
    allAll(pp) = mean(corrt);

    figure
    gscatter(stashLocations(:,plotcols(1)),stashLocations(:,plotcols(3)),stashLabels(:,1))
    xlabel(key.c.all(plotcols(1)))
    ylabel(key.c.all(plotcols(3)))
    str = sprintf('leave segment out, donor %i',pp)
    title(str)
end


corrall = classStore(:,1) == classStore(:,2);
accall = mean(corrall)

corrgrasp = GraspClassStore(:,1) == GraspClassStore(:,2);
accgrasp = mean(corrgrasp)


figure
gscatter3(dataStore(:,plotcols(1)),dataStore(:,plotcols(2)),dataStore(:,plotcols(3)),classStore(:,1))
xlabel(key.c.all(plotcols(1)))
ylabel(key.c.all(plotcols(2)))
zlabel(key.c.all(plotcols(3)))
title('Est Class Intra')

