%test window discriminant
clear all



%% 2D

nn = 1000;
az = 20;
el = 30;


MU1 = [0 0];
SIGMA1 = [3 0.1; 0.1 0.3];

MU2 = [0 0];
SIGMA2 = [0.3 0; 0 3];

MU3 = [5 0];
SIGMA3 = [2 0; 0 1];

X1 = mvnrnd(MU1,SIGMA1,nn);
X2 = mvnrnd(MU2,SIGMA2,nn);
%X3 = mvnrnd(MU3,SIGMA3,nn);

Data = [X1;X2];
Y = [ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];
classes = unique(Y);

figure
scatter(X1(:,1),X1(:,2),'r.')
hold on
scatter(X2(:,1),X2(:,2),'b.')
%hold on
%scatter(X3(:,1),X3(:,2),'g.')
hold off
title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast');
set(hl,'FontSize',12);

filename = 'C:\Users\MRDLAB\Documents\Rod Dockter\Surgical Robotics\Dockter_PHD\Thesis\Figures\gaussdatasimple'
%PrettyPlots(filename)

%% plot some shiz
h = windowDiscriminantPlot2(Data, Y);

filename = 'C:\Users\MRDLAB\Documents\Rod Dockter\Surgical Robotics\Dockter_PHD\Thesis\Figures\discrimwindow1'
%PrettyPlots(filename)


%% Train that shiz
pthresh = 0.3;
kmeanlimit = 4;

if(~exist('kmg'))
    kmg = [];
end
[Models, kmg] = windowDiscriminantTrain(Data, Y, pthresh, kmg, kmeanlimit, 'ploton');

filename = 'C:\Users\MRDLAB\Documents\Rod Dockter\Surgical Robotics\Dockter_PHD\Thesis\Figures\discrimwindow3'
%PrettyPlots(filename)

%% Classify with bulk data set

classthresh = pthresh;

Classify = [];

for kk = 1:length(classes)
    %Get seperate data for each class
    datac{kk} = Data(Y == classes(kk),:);
end

for ii = 1:length(Data)
    datatemp = Data(ii,:);
    Pt = [];
    for cc = 1:length(classes)
        for mm = 1:length(Models{cc}.cluster)
            %%get scaled probability for each model at each data point
            Pt(mm,cc) = Models{cc}.cluster{mm}.scale * gaussianProbMV(datatemp,Models{cc}.cluster{mm}.sigma,Models{cc}.cluster{mm}.mu);
        end
    end
    valt(ii,1) = max(Pt(:));
    [rowm,colm] = find(Pt == valt(ii));
    if(valt(ii)) >= classthresh
        Classify(ii,:) = [Y(ii),colm];
    else
        Classify(ii,:) = [Y(ii),0];
    end
end


figure
scatter(datac{1}(:,1),datac{1}(:,2),'r.')
hold on
scatter(datac{2}(:,1),datac{2}(:,2),'b.')
hold on
handle = Surface3D(Data(:,1),Data(:,2),valt);
hold off
title('Classification Probability')
xlabel('x1','FontSize', 12)
ylabel('x2','FontSize', 12)
zlabel('Probability','FontSize', 12)
hl = legend('class1','class2')
set(hl,'FontSize',12);

filename = 'C:\Users\MRDLAB\Documents\Rod Dockter\Surgical Robotics\Dockter_PHD\Thesis\Figures\discrimwindowsep'
%PrettyPlots(filename)

idx = find(Classify(:,2) == 0);
temp = Classify;
temp(idx,:) = [];
Correct = temp(:,1) == temp(:,2);
Accuracy = sum(Correct)/length(Correct)

%% New data for testing

rng(42)

Xn1 = mvnrnd(MU1,SIGMA1,nn);
Xn2 = mvnrnd(MU2,SIGMA2,nn);
%X3 = mvnrnd(MU3,SIGMA3,nn);

DataNew = [Xn1;Xn2];
YNew = [ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];

figure
scatter(Xn1(:,1),Xn1(:,2),'r.')
hold on
scatter(Xn2(:,1),Xn2(:,2),'b.')
hold off
title('Online Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast')
set(hl,'FontSize',12);

%% Log Likelihood ratios

disp('classification accuracy with Rods Discriminant:')
LLR_X1 = windowDiscriminantOnline(Xn1, Models, classes, pthresh);
LLR_X1(end);
class_okay1 = LLR_X1 >= 0;
Accuracy1 = sum(class_okay1)/length(class_okay1)

LLR_X2 = windowDiscriminantOnline(Xn2, Models, classes, pthresh);
LLR_X2(end);
class_okay2 = LLR_X2 <= 0;
Accuracy2 = sum(class_okay2)/length(class_okay2)



%% Compare with vanilla gaussians



[NN,SS] = size(DataNew);

mu1 = mean(Xn1);
sig1 = cov(Xn1);
mu2 = mean(Xn2);
sig2 = cov(Xn2);

Pall = [];
for ii = 1:length(Xn1)
    dtemp = Xn1(ii,:); 
    ptg(1) = gaussianProbMV(dtemp,sig1,mu1);
    ptg(2) = gaussianProbMV(dtemp,sig2,mu2);
    
    Pall(ii,:) = ptg;
    
    if(ii > 2*SS) %once we have enough data
        [LLRG1(ii),combos] = LLRmulti(Pall);
    else
        LLRG1(ii) = 0; %we just dont know yet
    end
    
    pdiff1(ii,:) = ptg(1) - ptg(2);
    
    [m,idx] = max(ptg);
    
    classin1(ii,:) = [1, idx];
end

Pall = [];
for ii = 1:length(Xn2)
    dtemp = Xn2(ii,:); 
    ptg(1) = gaussianProbMV(dtemp,sig1,mu1);
    ptg(2) = gaussianProbMV(dtemp,sig2,mu2);
    
    Pall(ii,:) = ptg;
    
    if(ii > 2*SS) %once we have enough data
        [LLRG2(ii),combos] = LLRmulti(Pall);
    else
        LLRG2(ii) = 0; %we just dont know yet
    end
    
    pdiff2(ii,:) = ptg(1) - ptg(2);
    
    [m,idx] = max(ptg);
    
    classin2(ii,:) = [2, idx];
    
end


corr1lr = LLRG1 >= 0;
corr2lr = LLRG2 <= 0;

corr1 = classin1(:,1) == classin1(:,2);
corr2 = classin2(:,1) == classin2(:,2);

disp('classification accuracy with gaussians:')
accg1 = sum(corr1)/length(corr1)
accg2 = sum(corr2)/length(corr2)



%% Compare plots

figure
plot([1:length(LLR_X1)],LLR_X1,'r-*' )
hold on
plot([1:length(pdiff1)],pdiff1,'b--+' )
hold off
title('Discriminant Window vs. Gaussian Class 1','FontSize', 12)
xlabel('x1','FontSize', 12)
ylabel('x2','FontSize', 12)
hl = legend('LLR Classification','Gaussian Probabilities','Location','northeast');
set(hl,'FontSize',12);


figure
plot([1:length(LLR_X2)],LLR_X2,'g-*' )
hold on
plot([1:length(pdiff2)],pdiff2,'c--+' )
hold off
title('Discriminant Window vs. Gaussian Class 2','FontSize', 12)
xlabel('x1','FontSize', 12)
ylabel('x2','FontSize', 12)
hl = legend('LLR Classification','Gaussian Probabilities','Location','northeast');
set(hl,'FontSize',12);




