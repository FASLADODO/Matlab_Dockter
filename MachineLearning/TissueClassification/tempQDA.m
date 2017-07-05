%% load up da data

load tissueData.mat

%% (intra-patient) get training/validation data and classify

NG = length(PatientData{1}.tissue{1}.grasp);

plotvid = randi([1,NG],1);

looaccIntra = [];
labelMaker = [0;1];


%get training and validation data
for pp = 1:length(PatientData)
        
    for loo = 1:NG
        %preallocate
        DataTrain = [];
        LabelsTrain = [];
        
        %decide on our leave one out
        trainIDX = 1:NG;
        trainIDX(loo) = [];
        valIDX = loo;

        
        for gg = 1:length(PatientData{pp}.tissue)
            for ss = 1:length(PatientData{pp}.tissue{gg}.grasp)
                struct = PatientData{pp}.tissue{gg}.grasp{ss};
                
                if(any(trainIDX==ss))
                    DataTrain = [DataTrain; struct.State, struct.Input];
                    LabelsTrain = [LabelsTrain; ones(length(struct.Input),1)*labelMaker(gg) ];
                elseif(any(valIDX==ss))
                    DataValidate{gg} = [struct.State, struct.Input];
                    LabelsValidate{gg} = [ones(length(struct.Input),1)*labelMaker(gg) ];
                end
            end
        end
    

        %Now Try and get logistic regression params
        ModelLR =  LogRegTrain(DataTrain,LabelsTrain);

        %now validate
        for vv = 1:length(key.cslist)
            [P,L] = LogRegOnline(DataValidate{vv},ModelLR );

            %time stamp accuracy
            classtemp = L > 0.5;
            corrtemp = LabelsValidate{vv} == classtemp;
            acctemp = mean(corrtemp);
            looaccIntra = [looaccIntra; acctemp];

        end

        %for plotsz
        if(any(valIDX==plotvid))
            stankd = [];
            for vv = 1:length(key.cslist)
                [P,L] = LogRegOnline(DataValidate{vv},ModelLR );
                stankd =[stankd; DataValidate{vv}, L];
            end
            figure
            scatter3(stankd(:,key.theta),stankd(:,key.thetadot),stankd(:,key.thetadotdot),10,stankd(:,end));
            xlabel('theta')
            ylabel('thetadot')
            zlabel('thetadotdot')
            colorbar
            colormap winter
            title('logistic regression mapping')
        end
    end
end


disp('intra patient accuracy')
mean(looaccIntra)

%% (inter-patient) get training/validation data and classify

NP = length(PatientData);

plotvid = randi([1,NP],1);

looaccInter = [];
labelMaker = [0;1];


%get training and validation data
for loo = 1:NP
    
    %preallocate
    DataTrain = [];
    LabelsTrain = [];

    %decide on our leave one out
    trainIDX = 1:NP;
    trainIDX(loo) = [];
    valIDX = loo;

    
    for pp = 1:length(PatientData)
        for gg = 1:length(PatientData{pp}.tissue)
            for ss = 1:length(PatientData{pp}.tissue{gg}.grasp)
                struct = PatientData{pp}.tissue{gg}.grasp{ss};
                
                if(any(trainIDX==pp))
                    DataTrain = [DataTrain; struct.State, struct.Input];
                    LabelsTrain = [LabelsTrain; ones(length(struct.Input),1)*labelMaker(gg) ];
                elseif(any(valIDX==pp))
                    DataValidate{gg} = [struct.State, struct.Input];
                    LabelsValidate{gg} = [ones(length(struct.Input),1)*labelMaker(gg) ];
                end
            end
        end
    end
    

    %Now Try and get logistic regression params
    ModelLR =  LogRegTrain(DataTrain,LabelsTrain);

    %now validate
    for vv = 1:length(key.cslist)
        [P,L] = LogRegOnline(DataValidate{vv},ModelLR );

        %time stamp accuracy
        classtemp = L > 0.5;
        corrtemp = LabelsValidate{vv} == classtemp;
        acctemp = mean(corrtemp);
        looaccInter = [looaccInter; acctemp];

    end

    %for plotsz
    if(any(valIDX==plotvid))
        
        [P,L] = LogRegOnline(DataTrain,ModelLR );
        stankd =[DataTrain, L];
            
        figure
        scatter3(stankd(:,key.theta),stankd(:,key.thetadot),stankd(:,key.thetadotdot),10,stankd(:,end));
        xlabel('theta')
        ylabel('thetadot')
        zlabel('thetadotdot')
        colorbar
        colormap winter
        title('logistic regression mapping')
    end
end


disp('inter patient accuracy')
mean(looaccInter)


 %% create a logistic regression classifier

ModelLRAll =  LogRegTrain(AllData,AllLabels);

%% try classifying all data

[P,Lall] = LogRegOnline(AllData,ModelLRAll );

figure
scatter3(AllData(:,key.theta),AllData(:,key.thetadot),AllData(:,key.thetadotdot),10,Lall);
xlabel('theta')
ylabel('thetadot')
zlabel('thetadotdot')
colorbar
colormap winter
title('logistic regression mapping')

