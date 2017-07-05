%% load up da data

load segData.mat

Linear_Func = @(X) [X(:,key.thetadotdot), X(:,key.thetadot), X(:,key.theta), X(:,key.theta).^4];

cslist = {'liver'; 'pancreas'};

%% get true linear params

%simulink things
tend = 1;
T = 0.001; % sampling period is fronm 1KHz
t = 0:T:tend;
F = 4;
%input is linear increase with time
input.time = t;
input.signals.values = F*t; %10*ones(1,length(t));
input.time = [input.time]';
input.signals.values = [input.signals.values]';
input.signals.dimensions = 1;


for gg = 1:length(cslist)
    ClassData{gg} = [];
    ClassInput{gg} = [];
end
% get all data from all grasps as a test
for pp = 1:length(SegData.Donor)
    for gg = 1:length(SegData.Donor{pp}.Tissue)
        for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
            for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                struct = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss};

                ClassData{gg} = [ClassData{gg}; struct.State];
                ClassInput{gg} = [ ClassInput{gg}; struct.Input];
            end
        end
    end
end

%run with true params for all data
for gg = 1:length(cslist)
    dtemp = Linear_Func(ClassData{gg});
    utemp = ClassInput{gg};
    phi = pinv(dtemp)*utemp;
    % SIMULINK
    sim('TissueModel1_linear.slx');

    % no noise, get output stuff
    True_U{gg} = input_out.Data(:,1);
    True_State{gg} = state.Data;
end

figure
gscatter3(AllData(:,key.theta),AllData(:,key.thetadot),AllData(:,key.thetadotdot),AllLabels)
hold on
gg = 1;
scatter3(True_State{gg}(:,key.theta),True_State{gg}(:,key.thetadot),True_State{gg}(:,key.thetadotdot),'c.')
hold on
gg = 2;
scatter3(True_State{gg}(:,key.theta),True_State{gg}(:,key.thetadot),True_State{gg}(:,key.thetadotdot),'k.')
hold off
xlabel('theta')
ylabel('thetadot')


%% least squares classification with leave patient out (inter patient)

epochs = length(SegData.Donor);
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
    for gg = 1:length(cslist)
        DataTrain{gg} = [];
        InputTrain{gg} = [];
        LabelsTrain{gg} = [];
    end
    DataValidate = [];
    InputValidate = [];
    LabelsValidate = [];
    
    
    %get training and validation data
    for pp = 1:length(SegData.Donor)
        for gg = 1:length(SegData.Donor{pp}.Tissue)
            segid = 1;
            for ll = 1:length(SegData.Donor{pp}.Tissue{gg}.Location )
                for ss = 1:length(SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp )
                    struct = SegData.Donor{pp}.Tissue{gg}.Location{ll}.Grasp{ss};
                    
                    if(any(trainIDX==pp))
                        DataTrain{gg} = [DataTrain{gg}; struct.State];
                        InputTrain{gg} = [ InputTrain{gg}; struct.Input];
                    elseif(any(valIDX==pp))
                        DataValidate.Tissue{gg}.Grasp{segid} = struct.State;
                        InputValidate.Tissue{gg}.Grasp{segid} = struct.Input;
                        segid = segid + 1;
                    end
                end
            end
        end
    end

    %Train params with training set
    for gg = 1:length(key.cslist)
        dtemp =  Linear_Func(DataTrain{gg});
        utemp = InputTrain{gg};
        ParamsTrain{gg} = pinv(dtemp)*utemp;
    end
    
    %test classification with validation set
    for gg = 1:length(DataValidate.Tissue)
        for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
            dtemp = Linear_Func(DataValidate.Tissue{gg}.Grasp{ss});
            utemp = InputValidate.Tissue{gg}.Grasp{ss};
            error = [];
            for cc = 1:length(key.cslist)
                error = [error, abs(utemp - dtemp*ParamsTrain{cc})];
            end
            errorStore = [errorStore; error];
            [~,classest] = min(error,[],2);
            totalest = mean(classest);
            graspclass = (totalest > 1.5) + 1;
            GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
            classStore = [ classStore; classest, ones(length(error),1)*gg];
            dataStore = [dataStore; DataValidate.Tissue{gg}.Grasp{ss}, utemp];
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
for pp = 1:length(SegData.Donor)
    
    %decide on our leave one out based on largest number of locations
    for gg = 1:length(SegData.Donor{pp}.Tissue)
        epochall(gg) = length(SegData.Donor{pp}.Tissue{gg}.Location );
    end
    epochs=max(epochall);
        
    for oo = 1:epochs
        %for storage
        for gg = 1:length(cslist)
            DataTrain{gg} = [];
            InputTrain{gg} = [];
            LabelsTrain{gg} = [];
        end
        DataValidate = [];
        InputValidate = [];
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
                        DataTrain{gg} = [DataTrain{gg}; struct.State];
                        InputTrain{gg} = [ InputTrain{gg}; struct.Input];
                    elseif(any(valIDX==ll))
                        DataValidate.Tissue{gg}.Grasp{segid} = struct.State;
                        InputValidate.Tissue{gg}.Grasp{segid} = struct.Input;
                        segid = segid + 1;
                    end
                end
            end
        end


        %Train params with training set
        for gg = 1:length(cslist)
            dtemp =  Linear_Func(DataTrain{gg});
            utemp = InputTrain{gg};
            %least squares parameters
            ParamsTrain{gg} = pinv(dtemp)*utemp;
        end

        %test classification with validation set
        for gg = 1:length(DataValidate.Tissue)
            for ss = 1:length(DataValidate.Tissue{gg}.Grasp)
                dtemp = Linear_Func(DataValidate.Tissue{gg}.Grasp{ss});
                utemp = InputValidate.Tissue{gg}.Grasp{ss};
                %classify using least squares
                error = [];
                for cc = 1:length(key.cslist)
                    error = [error, abs(utemp - dtemp*ParamsTrain{cc})];
                end
                errorStore = [errorStore; error];
                [~,classest] = min(error,[],2);
                totalest = mean(classest);
                graspclass = (totalest > 1.5) + 1;
                GraspClassStore = [GraspClassStore; graspclass, gg, totalest];
                classStore = [ classStore; classest, ones(length(error),1)*gg];
                dataStore = [dataStore; DataValidate.Tissue{gg}.Grasp{ss}, utemp];
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
