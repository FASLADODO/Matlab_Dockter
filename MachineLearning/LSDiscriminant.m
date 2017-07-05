
runs = 50;

nn = 100;
noiz = 0.2;

rangen = [0 6*pi];

paramsg1 = [2; 3; -0.3; -0.8]; %m, b
paramsg2 = [2; 3; -0.2; -0.4]; %m, b

X1_train = [];
X2_train = [];

LSD_train = [];
LSY_train = [];
LABELS_train = [];


DataOn = [];

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
    
    X1_train = [X1_train; X1];
    X2_train = [X2_train; X2];

    LSD_train = [LSD_train; d1; d2];
    LSY_train = [LSY_train; yx1; yx2];
    LABELS_train = [LABELS_train; ones(nn,1)*0; ones(nn,1)*1];
    
    DataOn{1}{ii} = [d1, yx1];
    DataOn{2}{ii} = [d2, yx2];
end



figure
scatter(X1_train(:,1),X1_train(:,2),'r.')
hold on
scatter(X2_train(:,1),X2_train(:,2),'b.')
hold off
title('Sample Data','FontSize', 12)
xlabel('state1','FontSize', 12)
ylabel('state2','FontSize', 12)
hl = legend('class1','class2','Location','southeast');
set(hl,'FontSize',12);


%% Solve log reg and params


ModelFull = LLSRtrain(LSD_train,LSY_train,LABELS_train);


ip = 5;
%do it online
dtest = DataOn{1}{ip};

[P1,Lm1] = LogRegOnline(dtest,ModelFull.LogReg );
[LLR1,combos1] = cumSumLLR(P1);

figure
scatter(1:length(Lm1),Lm1,'b.')
title('logreg class 1')

figure
plot(P1)
title('Probabilities 1')

figure
scatter(1:length(LLR1),LLR1,'r.')
title('LLR class 1')

%do it online
dtest = DataOn{2}{ip};

[P2,Lm2] = LogRegOnline(dtest,ModelFull.LogReg  );
[LLR2,combos2] = cumSumLLR(P2);

figure
scatter(1:length(Lm2),Lm2,'b.')
title('logreg class 2')

figure
plot(P2)
title('Probabilities 2')

figure
scatter(1:length(LLR2),LLR2,'r.')
title('LLR class 2')

%% plot LS

classes = unique(LABELS_train);
np = 500;

spaceLS = [];

%get generic data matrix
xm = linspace(rangen(1),rangen(2),np)';
xm1 = xm;
spaceLS = [ones(np,1), xm1, xm1.^2, xm1.^3];

%map least squares params to get estimated trajectory
Y1 = spaceLS*ModelFull.LSParams{1};
Y2 = spaceLS*ModelFull.LSParams{2};

colordata = [];
%Get coloring
for cc = 1:length(classes)
    colorp{cc} = [];
    for ii = 1:length(DataOn{cc})

        colordata = DataOn{cc}{ii};

        [PALL,Lm] = LogRegOnline(colordata,ModelFull.LogReg  );
        [LLR,combos] = cumSumLLR(PALL);

        colorp{cc} = [colorp{cc};  abs(LLR)]; %(abs(Lm-0.5) * 2) ]; %PALL(:,cc)];
    end
end


%plot the fit
figure
scatter(X1_train(:,1),X1_train(:,2),2,colorp{1})
hold on
plot(spaceLS(:,2),Y1,'g-','LineWidth',2)
hold on
scatter(X2_train(:,1),X2_train(:,2),2,colorp{2})
hold on
plot(spaceLS(:,2),Y2,'c-','LineWidth',2)
hold off
colormap cool
colorbar

%% test seperation change

%find the diff in params to see seperation
differ = abs(Y1- Y2);
%fit the diff
p2 = polyfit(spaceLS(:,2),differ,3)

diffp = ModelFull.LSParams{2} - ModelFull.LSParams{1}

%plot the diff
figure
plot(spaceLS(:,2),differ,'r-')
hold on
plot(spaceLS(:,2),fliplr(spaceLS(:,[1,2,3,4]))*p2','b-')
hold off






