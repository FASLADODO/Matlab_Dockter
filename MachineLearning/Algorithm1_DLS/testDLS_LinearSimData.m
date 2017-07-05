% test classification using DLS on linear simulated data

nn = 1000;
xx = linspace(0,10,nn)'; 
Xall = [xx + randn(nn,1)*0.1, ones(nn,1)];

parm1 = [3;1]; %m, b
parm2 = [4;1];

Y1 = Xall * parm1 + randn(nn,1)*0.1;
Y2 = Xall * parm2 + randn(nn,1)*0.1;

X = [Xall; Xall];
Y = [Y1; Y2];
Labels = [ones(nn,1)*1; ones(nn,1)*2];

figure
gscatter(X(:,1),Y(:,1),Labels)


lmd = 0.4
tic
Parameters = DLS_Train(X, Y, Labels, 0.4)
toc

%% now lets classify

figure
for cc = 1:2
   xtemp = X(Labels == cc,:);
   ytemp = Y(Labels == cc,:);
   [ClassEst{cc}, ClassEstCum{cc}, Confidence{cc}] = DLS_Online(xtemp, ytemp, Parameters);

    
    gscatter(xtemp(:,1),ytemp,ClassEstCum{cc})
    hold on
end
hold off

figure
cc = 1;
plot(Confidence{cc}(:,1),'c')
hold on
plot(Confidence{cc}(:,2),'r')
hold off

%% load up linear data (start with kappa = 1.1)

load LinearData11.mat

%% grab all data for plots

coluse = [1,2,3];
ycol = 4;
mapping = [1,2];

XData = [];
YData = [];
Labels = [];


for cc = 1:length(SegmentData.Class)
    for ii = 1:length(SegmentData.Class{cc}.Iteration)
        [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});

        Labels = [Labels; ones(ns,1)*mapping(cc)];
        XData = [XData; SegmentData.Class{cc}.Iteration{ii}(:,coluse)];
        YData = [YData; SegmentData.Class{cc}.Iteration{ii}(:,ycol)];
    end
end


%% train time
tic
for t = 1:100
    params = DLS_Train(XData, YData, Labels, 0.1);
end
toc

%% online time
XDataOn = XData(1:1000,:);
YDataOn = YData(1:1000,:);
tic
for t = 1:100
    [cesttemp, cestcumtemp, conftemp, stemp] = DLS_Online(XDataOn, YDataOn, params);
end
toc

%% perform leave one out

fsize = 14;

% get data
classes = length(SegmentData.Class);
runs = length(SegmentData.Class{1}.Iteration );

%state and input columns
statecol = [1,2,3];
ycol = 4;

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
LambdaUsedAll = [];
ScoreAll = [];

%we also want to check which lambda is the best
LambdaAcc = [];
LambdaList = 0:0.1:1;

%leave one out epoch
for oo = 1:runs
    
    fprintf('Epoch %d of %d \n',oo,runs)
    
    %clear our matrices
    Xtrain = [];
    Ytrain = [];
    Labelstrain = [];
    Xtest = [];
    Ytest = [];
    Labelstest = [];
    Ptrain = []; %parameters

    %build up train and test data for this epoch
    for cc = 1:classes
        for ii = 1:runs
            %number of points in this iteration
            [ns,nc] = size(SegmentData.Class{cc}.Iteration{ii});
            
            if(ii == oo)
                %left out data
                Xtest = [Xtest; SegmentData.Class{cc}.Iteration{ii}(:,statecol)];
                Ytest = [Ytest; SegmentData.Class{cc}.Iteration{ii}(:,ycol)];
                Labelstest = [Labelstest; ones(ns,1)*cc ];
            else
                %training data
                Xtrain = [Xtrain; SegmentData.Class{cc}.Iteration{ii}(:,statecol)];
                Ytrain = [Ytrain; SegmentData.Class{cc}.Iteration{ii}(:,ycol)];
                Labelstrain = [Labelstrain; ones(ns,1)*cc ];
            end
        end
    end
    
    %figure out lambda
    %lambdatest = 0; %for now
    bestAcc = 0;
    ClassAllBest= [];
    ClassTimeBest = [];
    ConvTimeBest= [];
    ScoreBest = [];
    LambdaBest = 0;
    
    lid = 1;
    for lambdatest = LambdaList
    
        %now we train for this epoch
        Ptrain = DLS_Train(Xtrain, Ytrain, Labelstrain, lambdatest);

        %store classification
        ClassAllTemp = [];
        ClassTimeTemp = [];
        ConvTimeTemp = [];
        STimeTemp = [];
        
        %now we classify with test data
        for cc = 1:classes
            %grab one classes data
           xtemp = Xtest(Labelstest == cc,:);
           ytemp = Ytest(Labelstest == cc,:);
           labeltrue = ones(length(ytemp),1)*cc;

           %try classifying
           [cesttemp, cestcumtemp, conftemp, stemp] = DLS_Online(xtemp, ytemp, Ptrain);

           %find where the last time it changed classification was
           convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
           
           %store it temporarily
           ClassAllTemp = [ClassAllTemp; cc, cesttemp];
           ClassTimeTemp = [ClassTimeTemp; labeltrue, cestcumtemp];
           ConvTimeTemp = [ConvTimeTemp; convtime];
           STimeTemp = [STimeTemp; stemp];
        end
        
        %check if this lambda is best
        corrtemp = ClassTimeTemp(:,1) == ClassTimeTemp(:,2);
        acctemp = mean(corrtemp);
        %store it if its best
        if(acctemp > bestAcc)
           bestAcc = acctemp;
           ClassAllBest = ClassAllTemp;
           ClassTimeBest = ClassTimeTemp;
           ConvTimeBest = ConvTimeTemp;
           LambdaBest = lambdatest;
           ScoreBest =  STimeTemp;
        end
        
        %stash accuracy vs lambda for kicks
        LambdaAcc(oo,lid) = acctemp;
        lid = lid + 1;
    end
    
    %stash it
    ClassificationAll = [ClassificationAll; ClassAllBest ];
    ClassificationTimeAll = [ClassificationTimeAll; ClassTimeBest];
    ConvergenceTime = [ConvergenceTime; ConvTimeBest];
    LambdaUsedAll = [LambdaUsedAll; LambdaBest];
    ScoreAll = [ ScoreAll; ScoreBest];
end

avgconvtime = mean(ConvergenceTime)

corrtime = ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2);
acctime = mean(corrtime)
corr = ClassificationAll(:,1) == ClassificationAll(:,2);
acc = mean(corr)

%kappa=1.1 gets 92.9% time, 100%overall classification, 12 % convergence
%time
%kappa=1.2 gets 100% classification

meanacck11 = mean(LambdaAcc);
% meanacck12 = mean(LambdaAcc);

%% Scores for correct and incorrect classifications

incorrectScores = abs(ScoreAll(ClassificationTimeAll(:,1) ~= ClassificationTimeAll(:,2) ));
correctScores = abs(ScoreAll(ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2) ));
C = [ones(length(incorrectScores),1)*1; ones(length(correctScores),1)*2] ;

mean(incorrectScores)
std(incorrectScores)
mean(correctScores)
std(correctScores)

figure
boxplot([incorrectScores; correctScores],C,'notch','on','labels',{'Incorrect','Correct'})
xlabel('Classification','FontSize',fsize)
ylabel('Weight \alpha','FontSize',fsize)
set(findobj(gca,'Type','text'),'FontSize',fsize)

%% now we want to plot lambda vs accuracy for both kappas
%have to run previous cells twice and save accs for this to work

figure
plot(LambdaList,meanacck11*100,'c-o','LineWidth',3)
hold on
plot(LambdaList,meanacck12*100,'g-x','LineWidth',3)
hold off
xlabel('\lambda','FontSize',fsize)
ylabel('% Correct Classification','FontSize',fsize)
h_legend=legend('\kappa = 1.1', '\kappa = 1.2');
set(h_legend,'FontSize',fsize);
ylim([0,100])



save('lambdaaccslinear.mat','meanacck11','meanacck12')


%% For replotting

load 'lambdaaccslinear.mat'
LambdaList = 0:0.1:1;
fsize = 14;

figure
plot(LambdaList,meanacck11*100,'c-o','LineWidth',3)
hold on
plot(LambdaList,meanacck12*100,'g-x','LineWidth',3)
hold off
xlabel('\lambda','FontSize',fsize)
ylabel('% Correct Classification','FontSize',fsize)
h_legend=legend('\kappa = 1.1', '\kappa = 1.2');
set(h_legend,'FontSize',fsize);
ylim([0,100])






