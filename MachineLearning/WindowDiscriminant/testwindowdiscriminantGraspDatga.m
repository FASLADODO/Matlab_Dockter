clear all

%% fake grasper data runs
runs = 50;

nn = 100;
noiz = 0.2;

rangen = [0 6*pi];

paramsg1 = [2; 3; -0.3; -0.6]; %m, b
paramsg2 = [2; 3; -0.2; -0.4]; %m, b

Data_Online = [];
Y_Online = [];
Data_train = [];
Y_train = [];

X1_train = [];
X2_train = [];
X1_Online = [];
X2_Online = [];

LSD_online = [];
LSY_online = [];
LABELS_online = [];

LSD_train = [];
LSY_train = [];
LABELS_train = [];

for ii = 1:runs
    rng('shuffle')
    params1_temp = paramsg1 + randn(4,1)*0.01;
    params2_temp = paramsg2 + randn(4,1)*0.01;
    
    xm = linspace(rangen(1),rangen(2),nn)';

    xm1 = xm + randn(nn,1)*noiz;
    xm2 = xm + randn(nn,1)*noiz;

    d1 = [ones(nn,1), xm1, xm1.^2, xm1.^3];
    d2 = [ones(nn,1), xm2, xm2.^2, xm2.^3];

    yx1 = d1*params1_temp + randn(nn,1)*noiz;
    yx2 = d2*params2_temp + randn(nn,1)*noiz;


    X1 = [xm2, yx1];
    X2 = [xm1, yx2];
    

    
    if(ii < runs/2)
        Data_train = [Data_train; X1;X2];
        Y_train = [Y_train; ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];
        X1_train = [X1_train; X1];
        X2_train = [X2_train; X2];
        
        LSD_train = [LSD_train; d1;d2];
        LSY_train = [LSY_train; yx1; yx2];
        LABELS_train = [LABELS_train; ones(nn,1)*1;ones(nn,1)*2];
    else
        Data_Online = [Data_Online; X1;X2];
        Y_Online = [Y_Online; ones(nn,1)*1;ones(nn,1)*2]; %;ones(nn,1)*3];
        X1_Online = [X1_Online; X1];
        X2_Online = [X2_Online; X2];
        
        LSD_online = [LSD_online; d1;d2];
        LSY_online = [LSY_online; yx1; yx2];
        LABELS_online = [LABELS_online; ones(nn,1)*1;ones(nn,1)*2];
    end
    
end

classes = unique(Y_train);

figure
scatter(X1_train(:,1),X1_train(:,2),'r.')
hold on
scatter(X2_train(:,1),X2_train(:,2),'b.')
%hold on
%scatter(X3(:,1),X3(:,2),'g.')
hold off
title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast');
set(hl,'FontSize',12);


%% plot the resulting seperations
tic
[h] = windowDiscriminantPlot(Data_train, Y_train, 0.5);
toc

%% Train the models for each
pthresh = 0.5;
kmeanlimit = 7;

if(~exist('kmg'))
    kmg = [];
end
[Models, kmg] = windowDiscriminantTrain2(Data_train, Y_train, pthresh, kmg, kmeanlimit, 'ploton');



%% Log Likelihood ratios with online data

disp('classification accuracy with Rods Discriminant:')
LLR_X1 = windowDiscriminantOnline(X1_Online, Models, classes, pthresh);
LLR_X1(end);
class_okay1 = LLR_X1 >= 0;
Accuracy1 = sum(class_okay1)/length(class_okay1)

LLR_X2 = windowDiscriminantOnline(X2_Online, Models, classes, pthresh);
LLR_X2(end);
class_okay2 = LLR_X2 <= 0;
Accuracy2 = sum(class_okay2)/length(class_okay2)


params = LeastSqauresTrain(LSD_train,LSY_train,LABELS_train);
disp('classification accuracy with least squares:')
AccLC = LeastSqauresClassify( LSD_online, LSY_online, LABELS_online,params);

figure
plot(LLR_X1,'r');
figure
plot(LLR_X2,'b');


