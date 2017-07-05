%% load up da data

load segData.mat

cslist = {'liver'; 'pancreas'};
labelMaker = [-1;1]; %becasue class values for logistic regression are 0 and 1
rounds = 50; %boosting rounds
trainmethod = {'AdaBoostM1';'Tree'};
% trainmethod = {'Subspace';'KNN'};

%% get overall adaboost tree

DataAll = [];
LabelsAll = [];

for pp = 1:length(SegData.Donor)
    for gg = 1:length(SegData.Donor{pp}.Tissue)
        for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
            for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                struct = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss};

                DataAll = [DataAll; struct.State, struct.Input];
                LabelsAll = [LabelsAll; ones(length(struct.Input),1)*labelMaker(gg) ];
            end
        end
    end
end
    
%train adaboost
%TreeAll = AdaBoostTrain(DataAll,LabelsAll,rounds,false);
%now fit the adaboost ensemble tree
tic
TreeAll = fitensemble(DataAll,LabelsAll,trainmethod{1},rounds,trainmethod{2});
toc

%%
%estimate labels using the Tree
%[Est,Value] = AdaBoostClassify(DataAll,TreeAll);
tic
Est = predict(TreeAll,DataAll); 
LabelsEst = Est;
toc

%%

figure
gscatter3(DataAll(:,key.theta),DataAll(:,key.thetadot),DataAll(:,key.thetadotdot),LabelsAll)
xlabel('theta')
ylabel('thetadot')
title('true class')

% figure
% scatter3(DataAll(:,key.theta),DataAll(:,key.thetadot),DataAll(:,key.thetadotdot),15,Value)
% xlabel('theta')
% ylabel('thetadot')
% title('log output')
% colorbar
% colormap winter

figure
gscatter3(DataAll(:,key.theta),DataAll(:,key.thetadot),DataAll(:,key.thetadotdot),LabelsEst)
xlabel('theta')
ylabel('thetadot')
title('Est class')

corr = LabelsEst == LabelsAll;
acc = mean(corr)

%% adaboost classification with leave patient out (inter patient)

epochs = length(SegData.Donor);
trainIDX = 1:epochs;
valIDX = 1;

%storage
classStore = [];
GraspClassStore = [];
dataStore = [];

tic
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
                        LabelsTrain = [LabelsTrain; ones(length(struct.Input),1)*labelMaker(gg)  ];
                    elseif(any(valIDX==pp))
                        DataValidate.Tissue{gg}.Grasp{segid} = [struct.State, struct.Input];
                        LabelsValidate.Tissue{gg}.Grasp{segid} = [ones(length(struct.Input),1)*labelMaker(gg)  ];
                        segid = segid + 1;
                    end
                end
            end
        end
    end

    %Train params with training set
    clear TreeModel
    TreeModel = fitensemble(DataTrain,LabelsTrain,trainmethod{1},rounds,trainmethod{2});
    
    %test classification with validation set
    for gg = 1:length(DataValidate.Tissue)
        for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
            dtemp = DataValidate.Tissue{gg}.Grasp{ss};
            %predict using adaboost tree
            Est = predict(TreeModel,dtemp); 
            classest = Est;
            %overall estimate
            totalest = mean(Est);
            graspclass = sign(totalest);
            GraspClassStore = [GraspClassStore; graspclass, labelMaker(gg), totalest];
            classStore = [ classStore; classest, ones(length(classest),1)*labelMaker(gg) ];
            dataStore = [dataStore; dtemp];
        end
    end
end
toc

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

%%  adaboost  classification with leave location out (intra patient)

epochs = length(SegData.Donor{1}.Tissue{1}.Location);
trainIDX = 1:epochs;
valIDX = 1;

%storage
errorStore = [];
classStore = [];
GraspClassStore = [];
dataStore = [];

tic
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
                        LabelsTrain = [LabelsTrain; ones(length(struct.Input),1)*labelMaker(gg)  ];
                    elseif(any(valIDX==ll))
                        DataValidate.Tissue{gg}.Grasp{segid} = [struct.State, struct.Input];
                        LabelsValidate.Tissue{gg}.Grasp{segid} = [ones(length(struct.Input),1)*labelMaker(gg)  ];
                        segid = segid + 1;
                    end
                end
            end
        end


        %Train params with training set
        clear TreeModel
        TreeModel = fitensemble(DataTrain,LabelsTrain,trainmethod{1},rounds,trainmethod{2});



        %test classification with validation set
        for gg = 1:length(DataValidate.Tissue)
            for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
                dtemp = DataValidate.Tissue{gg}.Grasp{ss};
                %predict using adaboost tree
                Est = predict(TreeModel,dtemp); 
                classest = Est;
                %overall estimate
                totalest = mean(Est);
                graspclass = sign(totalest);
                GraspClassStore = [GraspClassStore; graspclass, labelMaker(gg), totalest];
                classStore = [ classStore; classest, ones(length(classest),1)*labelMaker(gg) ];
                dataStore = [dataStore; dtemp];
            end
        end       
    end
end
toc

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