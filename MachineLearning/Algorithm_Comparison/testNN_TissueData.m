% test classification using NN on tissue data

%load it up
load SmartToolSegments.mat

fsize = 14;


%% grab all data for plots

fsize = 14;

DataAll = [];
Labels = [];
mapping = [0,1];

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



%% perform leave one out

fsize = 14;

% get data
classes = 2;
runs = 5;
hiddenLayerSize = 5; %use 5 layers
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. NFTOOL falls back to this in low memory situations.
trainFcn = 'trainlm';  % Bayesian Regularization

%state and input columns
coluse = [key.c.Stress, key.c.Strain, key.c.dStrain];
mapping = [0,1];

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
ScoreTime = [];
DataStore = [];

%leave one out epoch
for oo = 1:runs
    
    fprintf('Epoch %d of %d \n',oo,runs)
    
    %clear our matrices
    Xtrain = [];
    Labelstrain = [];
    Xtest = [];
    Labelstest = [];
    
    %build up train and test data for this epoch
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
    
    %setup NN
    net = patternnet(hiddenLayerSize,trainFcn);

    net.divideParam.trainRatio = 75/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 10/100;
    net.trainParam.showWindow = false;
    net.trainParam.show=50;  %# of ephocs in display
    net.trainParam.lr=0.05;  %learning rate
    net.trainParam.epochs=100;  %max epochs
    net.trainParam.goal=0.05^2;  %training goal
    net.performFcn='mse';  %Name of a network performance function %type help nnperformance
    
    % Train the Network
    [net,~] = train(net,Xtrain',Labelstrain');

    %now we classify with test data
    for cc = 1:classes
        %grab one classes data
       xtemp = Xtest(Labelstest == mapping(cc),:);
       labeltrue = ones(size(xtemp,1),1)*mapping(cc);

       %try classifying 
       esttime = net(xtemp');
       esttime = esttime'; %transpose
       
       %over all
       cesttemp = (mean(esttime) >= 0.5);
       
       
       %cumilative
       cestcumtemp = (CumMean(esttime) >= 0.5);
       
       %score
       scoretimetemp = esttime;
        
        
       %find where the last time it changed classification was
       convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
       if(isempty(convtime))
           convtime = 0;
       end
       
       %store it temporarily
       ClassificationAll = [ClassificationAll; mapping(cc), cesttemp];
       ClassificationTimeAll = [ClassificationTimeAll; labeltrue, cestcumtemp];
       ConvergenceTime = [ConvergenceTime; convtime, mapping(cc)];
       DataStore = [DataStore; xtemp];
       ScoreTime.Class{cc}.run{oo} = scoretimetemp;
    end
        
    clear net
end

%convergence time
avgconvtime = mean(ConvergenceTime(:,1))
idxconv1 = find(ConvergenceTime(:,2) == 0);
idxconv2 = find(ConvergenceTime(:,2) == 1);
avgconvtime1 = mean(ConvergenceTime(idxconv1,1))
avgconvtime2 = mean(ConvergenceTime(idxconv2,1))

%inidividual accs (1: liver, 2 : pancreas)
idxTime1 = find(ClassificationTimeAll(:,1) == 0);
idxTime2 = find(ClassificationTimeAll(:,1) == 1);
corrtime1 = ClassificationTimeAll(idxTime1,1) == ClassificationTimeAll(idxTime1,2);
corrtime2 = ClassificationTimeAll(idxTime2,1) == ClassificationTimeAll(idxTime2,2);
acctime1 = mean(corrtime1)
acctime2 = mean(corrtime2)
idx1 = find(ClassificationAll(:,1) == 0);
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

%% leave one Location out (intra-patient variability)

% get data
classes = 2;
runs = 5;


%state and input functions (xdot + x + x^2
coluse = [key.c.Stress, key.c.Strain, key.c.dStrain];
mapping = [0,1];

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
ScoreTime = [];
DataStore = [];


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

        %setup NN
        net = patternnet(hiddenLayerSize,trainFcn);

        net.divideParam.trainRatio = 75/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 10/100;
        net.trainParam.showWindow = false;
        net.trainParam.show=50;  %# of ephocs in display
        net.trainParam.lr=0.05;  %learning rate
        net.trainParam.epochs=100;  %max epochs
        net.trainParam.goal=0.05^2;  %training goal
        net.performFcn='mse';  %Name of a network performance function %type help nnperformance

        % Train the Network
        [net,~] = train(net,Xtrain',Labelstrain');


        %now we classify with test data
        for cc = 1:classes
           %grab one classes data
           xtemp = Xtest(Labelstest ==  mapping(cc),:);
           labeltrue = ones(size(xtemp,1),1)*mapping(cc);

           %try classifying 
           esttime = net(xtemp');
           esttime = esttime'; %transpose

           %over all
           cesttemp = (mean(esttime) >= 0.5);


           %cumilative
           cestcumtemp = (CumMean(esttime) >= 0.5);

           %score
           scoretimetemp = esttime;

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
        end
        clear net
    end
end

%convergence time
avgconvtime = mean(ConvergenceTime(:,1))
idxconv1 = find(ConvergenceTime(:,2) == 0);
idxconv2 = find(ConvergenceTime(:,2) == 1);
avgconvtime1 = mean(ConvergenceTime(idxconv1,1))
avgconvtime2 = mean(ConvergenceTime(idxconv2,1))

%inidividual accs (1: liver, 2 : pancreas)
idxTime1 = find(ClassificationTimeAll(:,1) == 0);
idxTime2 = find(ClassificationTimeAll(:,1) == 1);
corrtime1 = ClassificationTimeAll(idxTime1,1) == ClassificationTimeAll(idxTime1,2);
corrtime2 = ClassificationTimeAll(idxTime2,1) == ClassificationTimeAll(idxTime2,2);
acctime1 = mean(corrtime1)
acctime2 = mean(corrtime2)
idx1 = find(ClassificationAll(:,1) == 0);
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

