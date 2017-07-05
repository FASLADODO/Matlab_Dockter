% test classification using NN on linear simulated data

% load up linear data

load LinearData11.mat

fsize = 14;

%% grab all data for plots

coluse = [1,4];
mapping = [0,1];

Data = [];
Labels = [];


for cc = 1:length(SegmentData.Class)
    for ii = 1:length(SegmentData.Class{cc}.Iteration)
        [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});
        
        Data = [Data; SegmentData.Class{cc}.Iteration{ii}(:,coluse)];
        Labels = [Labels; ones(ns,1)*mapping(cc)];
    end
end

figure
gscatter(Data(:,1),Data(:,2),Labels,'rc')
title('true class')



%% setup network

hiddenLayerSize = 5; %use 5 layers
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. NFTOOL falls back to this in low memory situations.
trainFcn = 'trainlm';  % Bayesian Regularization

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

%% train time

tic
for t = 1:10
    t
    % Train the Network
    [net,~] = train(net,Data',Labels');
end
toc

%% online time
DataOn = Data(1:1000,:);
tic
for t = 1:10
    t
    esttime = net(DataOn');
end
toc

%% Train the Network

[net,~] = train(net,Data',Labels');


outputs = net(Data');

figure
scatter(Data(:,1),Data(:,2),10,outputs)
title('est class')
colorbar
colormap cool

%% perform leave one out

fsize = 14;

% get data
classes = length(SegmentData.Class);
runs = length(SegmentData.Class{1}.Iteration );

%state and input columns
statecol = [1,2,4];
mapping = [0,1];

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
ScoreTime = [];
DataAll = [];

%leave one out epoch
for oo = 1:runs
    
    fprintf('Epoch %d of %d \n',oo,runs)
    
    %clear our matrices
    Xtrain = [];
    Labelstrain = [];
    Xtest = [];
    Labelstest = [];
    
    %build up train and test data for this epoch
    for cc = 1:classes
        for ii = 1:runs
            [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});
        
            
            %its too damn big
            ssmp = 1:100:ns;
            
            tempstate = SegmentData.Class{cc}.Iteration{ii}(ssmp,statecol);
            
            
            if(ii == oo)
                %left out data
                Xtest = [Xtest; tempstate];
                Labelstest = [Labelstest; ones(length(ssmp),1)*mapping(cc) ];
            else
                %training data
                Xtrain = [Xtrain; tempstate];
                Labelstrain = [Labelstrain; ones(length(ssmp),1)*mapping(cc) ];
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
       ConvergenceTime = [ConvergenceTime; convtime];
       DataAll = [DataAll; xtemp];
       ScoreTime.Class{cc}.run{oo} = scoretimetemp;
    end
        
    clear net
end

avgconvtime = mean(ConvergenceTime)

corrtime = ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2);
acctime = mean(corrtime)
corr = ClassificationAll(:,1) == ClassificationAll(:,2);
acc = mean(corr)

%kappa=1.1 gets 72% time, 100%overall, 28 % convergence
%time
%kappa=1.2 gets 75% time, 100% classification, 36% convergence


%% plot confidence over time too

figure
for cc = 1:classes
    for oo = 1:runs
        stemp = ScoreTime.Class{cc}.run{oo};
        nns = length(stemp);
        xper = ([1:nns] ./ nns )*100;
        if(mapping(cc) == -1)
            h1=plot(xper,stemp,'r-');
            hold on
        else
            h2=plot(xper,stemp,'c-');
            hold on
        end
    end
end
hold off
xlabel('% Trajectory','FontSize',fsize)
ylabel('L_{on}','FontSize',fsize)
lh=legend([h1(1),h2(1)],'-1','1');
lh.FontSize = fsize;
lh.Location = 'NorthWest';

