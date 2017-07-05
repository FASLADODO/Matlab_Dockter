%% load up da data

load tissueData.mat

NonLinear_Fun_1 = @(X,Param) Param(1)*X(:,1) + Param(2)*X(:,2) + Param(3)*exp(-Param(4)*X(:,3) );
tend = 1;

%% get true nonlinear params

paramInit{1} = [1,0.02,0.1,0.01];
paramInit{2} = [1,0.85,0.1,0.4];

Param_NL = [];

fprintf('Train Nonlinear params \n')
for gg = 1:length(key.cslist)
    dtemp = ClassData{gg};
    utemp = ClassInput{gg};
    
    %train the paraemeters using gradient descent
    paramguess = paramInit{gg}';
    paramrefine = gradientdescentNL(NonLinear_Fun_1, paramguess, utemp, dtemp, 0.2, 200, 'batch', [1,2,3,3]);
    paramrefine = paramrefine./paramrefine(1);
    
    %change to be [d,a,b]
    phi = paramrefine([2,3,4])
    
    % SIMULINK
    sim('TissueModel1.slx');

    % no noise, get output stuff
    NL_U{gg} = input_out.Data(:,1);
    NL_State{gg} = state.Data;
    
    % GEt standard deviations in linear gaussian
    uest = NonLinear_Fun_1(dtemp,paramrefine);
    restemp = abs(utemp - uest);
    
    %store these things
    Param_NL{gg}.params = paramrefine;
    Param_NL{gg}.sigma = std(restemp);
    
end

figure
gscatter(AllData(:,key.theta),AllData(:,key.thetadot),AllLabels)
hold on
gg = 1;
scatter(NL_State{gg}(:,key.theta),NL_State{gg}(:,key.thetadot),'g.')
hold on
gg = 2;
scatter(NL_State{gg}(:,key.theta),NL_State{gg}(:,key.thetadot),'k.')
hold off
xlabel('theta')
ylabel('thetadot')


%% Classify with non linear model
%storage

errorStore = [];
errorClass = [];
classStore = [];
dataStore = [];
    
fprintf('Classify w/ Nonlinear params \n')
%test classification with validation set
for gg = 1:length(key.cslist)
    dtemp = ClassData{gg};
    utemp = ClassInput{gg};
    
    error = [];
    %test classification with all models
    for cc = 1:length(key.cslist)
        uest = NonLinear_Fun_1(dtemp,Param_NL{cc}.params );
        error = [error, abs(utemp - uest)];
    end
    errorStore = [errorStore; error];
    errorClass{gg} = error(:,gg);
    
    %determine class from error
    [~,classest] = min(error,[],2);
    classStore = [ classStore; classest, ones(length(error),1)*gg];
    dataStore = [dataStore; dtemp, utemp];
end


corr = classStore(:,1) == classStore(:,2);
acc = mean(corr)

figure
gscatter(dataStore(:,key.theta),dataStore(:,key.thetadot),classStore(:,1))
xlabel('theta')
ylabel('thetadot')

figure
gg = 1;
scatter(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),10,errorClass{gg});
hold on
gg = 2;
scatter(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),10,errorClass{gg});
hold off
xlabel('theta')
ylabel('thetadot')
colorbar
colormap cool
title('error')

figure
gscatter3(dataStore(:,key.theta),dataStore(:,key.thetadot),dataStore(:,end),classStore(:,1))
xlabel('theta')
ylabel('thetadot')
zlabel('thetadotdot')

figure
gscatter(dataStore(:,key.theta),dataStore(:,end),classStore(:,1))
xlabel('theta')
ylabel('force')


%% non-linear gaussian
%storage

probStore = [];
diffStore = [];
dataStore = [];
diffStoreAll = [];
    
fprintf('Non-linear gaussian \n')
%test classification with validation set
for gg = 1:length(key.cslist)
    dtemp = ClassData{gg};
    utemp = ClassInput{gg};
    
    %test classification with all models
    for cc = 1:length(key.cslist)
        [PL,Uest] = NonLinearGaussian(dtemp,utemp,Param_NL{cc},NonLinear_Fun_1);
        if(gg == cc)
            pwithin = PL;
        else
            pbetween = PL;
        end
    end
    Diff = ComputeRBFDifference(pwithin,pbetween);

    %stash
    probStore{gg} = pwithin;
    diffStore{gg} = Diff;
    dataStore = [dataStore; dtemp];
    diffStoreAll = [diffStoreAll; Diff];
end


%figure out equation to map theta, thetadot to difference
mapFunc = @(X) [ X(:,key.theta) , X(:,key.thetadot)];
modelDiff = RelativeGuassianTrainDiffs(dataStore,diffStoreAll,mapFunc);

% linear spacing for theta, thetadot
nn = 100;
thetalin = linspace(min(dataStore(:,key.theta)),max(dataStore(:,key.theta)),nn)';
thetadotlin = linspace(min(dataStore(:,key.thetadot)),max(dataStore(:,key.thetadot)),nn)';
datalin = [thetalin, thetadotlin];
difflin = datalin*modelDiff.Params;


figure
gg = 1;
scatter3(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),diffStore{gg},'b.')
hold on
gg = 2;
scatter3(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),diffStore{gg},'r.')
hold on
scatter3(datalin(:,1),datalin(:,2),difflin,'k.')
hold off
xlabel('theta')
ylabel('thetadot')
zlabel('diff')


figure
gg = 1;
[~,diffscaled{gg}] = RelativeGuassianProcessOnline(ClassData{gg},modelDiff);
scatter3(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),diffStore{gg},20,diffscaled{gg})
hold on
gg = 2;
[~,diffscaled{gg}] = RelativeGuassianProcessOnline(ClassData{gg},modelDiff);
scatter3(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),diffStore{gg},20,diffscaled{gg})
hold on
scatter3(datalin(:,1),datalin(:,2),difflin,'k.')
hold off
xlabel('theta')
ylabel('thetadot')
zlabel('diff')
colorbar
colormap winter

figure
gg = 1;
scatter3(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),probStore{gg},'b.')
hold on
gg = 2;
scatter3(ClassData{gg}(:,key.theta),ClassData{gg}(:,key.thetadot),probStore{gg},'r.')
hold off
xlabel('theta')
ylabel('thetadot')
zlabel('probability')


%% try classifying with non linear relative gaussian

%now we use -1, 1 as our classes to make weighting easier
csLookup = [-1;1];

dataStore = [];
timeStore = [];
estStore = [];
classStore = [];

fprintf('Classify w/ Nonlinear Gaussian \n')
%test classification with validation set
for gg = 1:length(key.cslist)
    for ss = 1:length(GraspData{gg}.grasp)
        dtemp = GraspData{gg}.grasp{ss}.State;
        utemp = GraspData{gg}.grasp{ss}.Input;

        error = [];
        %test classification with all models
        for cc = 1:length(key.cslist)
            uest = NonLinear_Fun_1(dtemp,Param_NL{cc}.params );
            error = [error, abs(utemp - uest)];
        end

        %determine class from error
        [~,class12] = min(error,[],2);
        classEst = csLookup(class12);
        
        %determine weigthing from gaussian differences
        [~,diffWeight] = RelativeGuassianProcessOnline(dtemp,modelDiff);
        
        %weight our estimates from diff
        classWeight = classEst.*diffWeight;
        
        %stash all the things!
        dataStore = [dataStore; dtemp];
        timeStore = [ timeStore; GraspData{gg}.grasp{ss}.Time];
        estStore = [estStore; classWeight];
        classStore =[ classStore; ones(length(utemp),1)*gg];
    end
end


figure
gscatter(timeStore,estStore,classStore)
xlabel('theta')
ylabel('confidence')


