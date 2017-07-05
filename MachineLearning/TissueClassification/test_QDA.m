%% load up da data

load segData.mat

cslist = {'liver'; 'pancreas'};

%% get overall QDA params

DataAll = [];
LabelsAll = [];

for pp = 1:length(SegData.Donor)
    for gg = 1:length(SegData.Donor{pp}.Tissue)
        for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
            for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                struct = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss};

                DataAll = [DataAll; struct.State, struct.Input];
                LabelsAll = [LabelsAll; ones(length(struct.Input),1)*gg];
            end
        end
    end
end
    
%train QDA
MdlQuadLoo = fitcdiscr(DataAll,LabelsAll,'DiscrimType','quadratic');

%estimate labels
LabelsEst = predict(MdlQuadLoo,DataAll);

figure
gscatter3(DataAll(:,key.theta),DataAll(:,key.thetadot),DataAll(:,key.thetadotdot),LabelsAll)
xlabel('theta')
ylabel('thetadot')
title('true class')

figure
gscatter3(DataAll(:,key.theta),DataAll(:,key.thetadot),DataAll(:,key.thetadotdot),LabelsEst)
xlabel('theta')
ylabel('thetadot')
title('Est class')

%% QDA classification with leave patient out (inter patient)

epochs = length(SegData.Donor);
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
    for pp = 1:length(SegData.Donor)
        for gg = 1:length(SegData.Donor{pp}.Tissue)
            segid = 1;
            for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
                for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                    struct = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss};
                    
                    if(any(trainIDX==pp))
                        DataTrain = [DataTrain; struct.State, struct.Input];
                        LabelsTrain = [LabelsTrain; ones(length(struct.Input),1)*gg ];
                    elseif(any(valIDX==pp))
                        DataValidate.Tissue{gg}.Grasp{segid} = [struct.State, struct.Input];
                        LabelsValidate.Tissue{gg}.Grasp{segid} = [ones(length(struct.Input),1)*gg ];
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
            dtemp = DataValidate.Tissue{gg}.Grasp{ss};
            %predict using QDA
            classest = predict(ParamQuad,dtemp);
            %overall estimate
            totalest = mean(classest);
            graspclass = (totalest > 1.5) + 1;
            GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
            classStore = [ classStore; classest, ones(length(classest),1)*gg];
            dataStore = [dataStore; dtemp];
        end
    end
end


corrall = classStore(:,1) == classStore(:,2);
accall = mean(corrall)

corrgrasp = GraspClassStore(:,1) == GraspClassStore(:,2);
accgrasp = mean(corrgrasp)

figure
gscatter(dataStore(:,key.theta),dataStore(:,key.thetadot),classStore(:,1))
xlabel('theta')
ylabel('thetadot')
title('Est Class Inter')

figure
gscatter(dataStore(:,key.theta),dataStore(:,end),classStore(:,1))
xlabel('theta')
ylabel('force')
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
for pp = 1:length(SegData.Donor)
    
    %decide on our leave one out based on largest number of locations
    for gg = 1:length(SegData.Donor{pp}.Tissue)
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
        for gg = 1:length(SegData.Donor{pp}.Tissue)
            %the number of locations left out may change
            dropidx = mod(oo-1,epochall(gg))+1; %because some tissues may have fewer locations
            trainIDX = 1:epochall(gg);
            trainIDX(dropidx) = [];
            valIDX = dropidx;
            segid = 1;
            for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
                for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                    %grab the current data
                    struct = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss};

                    if(any(trainIDX==ll))
                        DataTrain = [DataTrain; struct.State, struct.Input];
                        LabelsTrain = [LabelsTrain; ones(length(struct.Input),1)*gg ];
                    elseif(any(valIDX==ll))
                        DataValidate.Tissue{gg}.Grasp{segid} = [struct.State, struct.Input];
                        LabelsValidate.Tissue{gg}.Grasp{segid} = [ones(length(struct.Input),1)*gg ];
                        segid = segid + 1;
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
                dtemp = DataValidate.Tissue{gg}.Grasp{ss};
                %predict using QDA
                classest = predict(ParamQuad,dtemp);
                %overall estimate
                totalest = mean(classest);
                graspclass = (totalest > 1.5) + 1;
                GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
                classStore = [ classStore; classest, ones(length(classest),1)*gg];
                dataStore = [dataStore; dtemp];
            end
        end       
    end
end


corrall = classStore(:,1) == classStore(:,2);
accall = mean(corrall)

corrgrasp = GraspClassStore(:,1) == GraspClassStore(:,2);
accgrasp = mean(corrgrasp)

figure
gscatter(dataStore(:,key.theta),dataStore(:,key.thetadot),classStore(:,1))
xlabel('theta')
ylabel('thetadot')
title('Est Class Intra')

figure
gscatter(dataStore(:,key.theta),dataStore(:,end),classStore(:,1))
xlabel('theta')
ylabel('force')
title('Est Class Intra')