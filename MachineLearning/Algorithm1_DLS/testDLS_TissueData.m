%Test DLS grid on smart tool data

%load it up
load SmartToolSegments.mat


%% run through all grasps and store

fsize = 14;

DataAll = [];
Labels = [];
map = [1,2];

%get grasp statistics

for ii = 1:5 %length(SegData.Donor)
    for jj = 1:2 %length(SegData.Donor{ii}.Tissue)
        for kk = 1:length(SegData.Donor{ii}.Tissue{jj}.Location)
            for ll = length(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp)
                nn = size(SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data , 1);
                DataAll = [DataAll; SegData.Donor{ii}.Tissue{jj}.Location{kk}.Grasp{ll}.Data ];
                Labels = [Labels; ones(nn,1)*map(jj) ];
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
gscatter3(DataAll(:,tc(1)), DataAll(:,tc(2)), DataAll(:,tc(3)),key.t.all(Labels),'rc')
xlabel(key.c.all(tc(1)),'FontSize',fsize)
ylabel(key.c.all(tc(2)),'FontSize',fsize)
[lh,ic,ip,it]=legend('show');
lh.FontSize = fsize;
lh.Location = 'NorthEast';


%% leave one patient out (inter patient variablitity)

% get data
classes = 2;
runs = 5;

%state and input functions (xdot + x + x^2
statefunc = @(X) [ X(:,key.c.dStrain), X(:,key.c.Strain), X(:,key.c.Strain).^2];
ycol = key.c.Stress;

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
LambdaUsedAll = [];
ScoreAll = [];


%we also want to check which lambda is the best
LambdaAcc = [];
LambdaList = 0:0.1:1;

for oo = 1:runs %length(SegData.Donor)
    
    fprintf('Epoch %d of %d \n',oo,runs)
    
    %clear our matrices
    Xtrain = [];
    Ytrain = [];
    Labelstrain = [];
    Xtest = [];
    Ytest = [];
    Labelstest = [];
    Ptrain = []; %parameters

    %stash just the training data
    for dd = 1:runs %length(SegData.Donor)
        for cc = 1:classes %tissue type
            for ll = 1:length(SegData.Donor{dd}.Tissue{cc}.Location) %locations
                for gg = length(SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp) %grasp
                    %get linearized version
                    tempstate = statefunc(SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp{gg}.Data);
                    tempu = SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp{gg}.Data(:,ycol);

                    %number of points in this iteration
                    [ns,nc] = size(tempstate);


                    if(dd == oo)
                        %left out data
                        Xtest = [Xtest; tempstate];
                        Ytest = [Ytest; tempu];
                        Labelstest = [Labelstest; ones(ns,1)*cc ];
                    else
                        %training data
                        Xtrain = [Xtrain; tempstate];
                        Ytrain = [Ytrain; tempu];
                        Labelstrain = [Labelstrain; ones(ns,1)*cc ];
                    end
                end
            end
        end
    end
    
    %figure out lambda
    %lambdatest = 0; %for now
    bestAcc = 0;
    ClassAllBest= [];
    ClassTimeBest = [];
    ConvTimeBest= [];
    LambdaBest = 0;
    ScoreBest = [];

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
           [cesttemp, cestcumtemp, conftemp,stemp] = DLS_Online(xtemp, ytemp, Ptrain);

           %find where the last time it changed classification was
           convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
           if(isempty(convtime) )
               convtime = 0;
           end
            
           %store it temporarily
           ClassAllTemp = [ClassAllTemp; cc, cesttemp];
           ClassTimeTemp = [ClassTimeTemp; labeltrue, cestcumtemp];
           ConvTimeTemp = [ConvTimeTemp; convtime, cc];
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
    
    figure
    gscatter(Xtrain(:,2),Ytrain(:,1),Labelstrain,'rc')
    hold on
    gscatter(Xtest(:,2),Ytest(:,1),Labelstest,'gk','ooo')
    hold off
    title('loo')

end

%convergence time
avgconvtime = mean(ConvergenceTime(:,1))
idxconv1 = find(ConvergenceTime(:,2) == 1);
idxconv2 = find(ConvergenceTime(:,2) == 2);
avgconvtime1 = mean(ConvergenceTime(idxconv1,1))
avgconvtime2 = mean(ConvergenceTime(idxconv2,1))

%inidividual accs (1: liver, 2 : pancreas)
idxTime1 = find(ClassificationTimeAll(:,1) == 1);
idxTime2 = find(ClassificationTimeAll(:,1) == 2);
corrtime1 = ClassificationTimeAll(idxTime1,1) == ClassificationTimeAll(idxTime1,2);
corrtime2 = ClassificationTimeAll(idxTime2,1) == ClassificationTimeAll(idxTime2,2);
acctime1 = mean(corrtime1)
acctime2 = mean(corrtime2)
idx1 = find(ClassificationAll(:,1) == 1);
idx2 = find(ClassificationAll(:,1) == 2);
corr1 = ClassificationAll(idx1,1) == ClassificationAll(idx1,2);
corr2 = ClassificationAll(idx2,1) == ClassificationAll(idx2,2);
acc1 = mean(corr1)
acc2 = mean(corr2)

%overall acc
corrtime = ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2);
acctime = mean(corrtime)
corr = ClassificationAll(:,1) == ClassificationAll(:,2);
acc = mean(corr)

lambdagood = mean(LambdaUsedAll)

%Leave one patient out gets

meanaccklopo = mean(LambdaAcc);


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


%% leave one Location out (intra-patient variability)

% get data
classes = 2;
runs = 5;

%state and input functions (xdot + x + x^2
statefunc = @(X) [ X(:,key.c.dStrain), X(:,key.c.Strain), X(:,key.c.Strain).^2];
ycol = key.c.Stress;

%store the classifications
ClassificationAll = [];
ClassificationTimeAll = [];
ConvergenceTime = [];
LambdaUsedAll = [];
ScoreAll = [];

%we also want to check which lambda is the best
LambdaAcc = [];
LambdaList = 0:0.1:1;

for dd = 1:runs %length(SegData.Donor)
    
    fprintf('Epoch %d of %d \n',dd,runs)
    
    %now we leave out one location per
    sublocation = max([length(SegData.Donor{dd}.Tissue{1}.Location),length(SegData.Donor{dd}.Tissue{2}.Location)]);
    for oo = 1:sublocation
    
        %clear our matrices
        Xtrain = [];
        Ytrain = [];
        Labelstrain = [];
        Xtest = [];
        Ytest = [];
        Labelstest = [];
        Ptrain = []; %parameters

        %stash just the training data
        for cc = 1:classes %tissue type
            %figure out which location we have to leave out
            %(different tissues have different # locations)
            numlocs = length(SegData.Donor{dd}.Tissue{cc}.Location);
            leaveout = mod(oo-1,numlocs)+1; %mod 1 indexed
            for ll = 1:numlocs %locations
                for gg = length(SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp) %grasp
                    %get linearized version
                    tempstate = statefunc(SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp{gg}.Data);
                    tempu = SegData.Donor{dd}.Tissue{cc}.Location{ll}.Grasp{gg}.Data(:,ycol);

                    %number of points in this iteration
                    [ns,nc] = size(tempstate);
                    
                    if(ll == leaveout)
                        %left out data
                        Xtest = [Xtest; tempstate];
                        Ytest = [Ytest; tempu];
                        Labelstest = [Labelstest; ones(ns,1)*cc ];
                    else
                        %training data
                        Xtrain = [Xtrain; tempstate];
                        Ytrain = [Ytrain; tempu];
                        Labelstrain = [Labelstrain; ones(ns,1)*cc ];
                    end
                end
            end
        end

        %figure out lambda
        %lambdatest = 0; %for now
        bestAcc = 0;
        ClassAllBest= [];
        ClassTimeBest = [];
        ConvTimeBest= [];
        LambdaBest = 0;
        ScoreBest = [];

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
               [cesttemp, cestcumtemp, conftemp,stemp] = DLS_Online(xtemp, ytemp, Ptrain);

               %find where the last time it changed classification was
               convtime = find(cestcumtemp ~= cestcumtemp(end),1,'last') /length(cestcumtemp);
               if(isempty(convtime) )
                   convtime = 0;
               end

               %store it temporarily
               ClassAllTemp = [ClassAllTemp; cc, cesttemp];
               ClassTimeTemp = [ClassTimeTemp; labeltrue, cestcumtemp];
               ConvTimeTemp = [ConvTimeTemp; convtime,cc];
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
end

%convergence time
avgconvtime = mean(ConvergenceTime(:,1))
idxconv1 = find(ConvergenceTime(:,2) == 1);
idxconv2 = find(ConvergenceTime(:,2) == 2);
avgconvtime1 = mean(ConvergenceTime(idxconv1,1))
avgconvtime2 = mean(ConvergenceTime(idxconv2,1))

%inidividual accs (1: liver, 2 : pancreas)
idxTime1 = find(ClassificationTimeAll(:,1) == 1);
idxTime2 = find(ClassificationTimeAll(:,1) == 2);
corrtime1 = ClassificationTimeAll(idxTime1,1) == ClassificationTimeAll(idxTime1,2);
corrtime2 = ClassificationTimeAll(idxTime2,1) == ClassificationTimeAll(idxTime2,2);
acctime1 = mean(corrtime1)
acctime2 = mean(corrtime2)
idx1 = find(ClassificationAll(:,1) == 1);
idx2 = find(ClassificationAll(:,1) == 2);
corr1 = ClassificationAll(idx1,1) == ClassificationAll(idx1,2);
corr2 = ClassificationAll(idx2,1) == ClassificationAll(idx2,2);
acc1 = mean(corr1)
acc2 = mean(corr2)

%overall acc
corrtime = ClassificationTimeAll(:,1) == ClassificationTimeAll(:,2);
acctime = mean(corrtime)
corr = ClassificationAll(:,1) == ClassificationAll(:,2);
acc = mean(corr)

lambdagood = mean(LambdaUsedAll)

%Leave one location out gets

meanaccklolo = mean(LambdaAcc);


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
plot(LambdaList,meanaccklolo*100,'c-o','LineWidth',3)
hold on
plot(LambdaList,meanaccklopo*100,'g-x','LineWidth',3)
hold off
xlabel('\lambda','FontSize',fsize)
ylabel('% Correct Classification','FontSize',fsize)
h_legend=legend('LOLO', 'LODO');
set(h_legend,'FontSize',fsize);
ylim([0,100])


save('lambdaaccstissue.mat','meanaccklolo','meanaccklopo')


