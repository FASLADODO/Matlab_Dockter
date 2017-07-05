%% THIS IS THE VERSION THAT USES FUNCTIONS

% USE IT TO TEST 1D, 2D, 3D and multi modal


%% Make some 1D rand data

dist = 3;

means{1} = [2];
sigmas{1} = [1];

means{2} = means{1} + dist;
sigmas{2} = [1];

nn = [5000,5000];

[Data,Labels] = CreateData(means,sigmas,nn);

figure
gscatter(Data(:,1),zeros(length(Data),1),Labels)

%Train the RBF model
[Fits] = RelativeRBFTrain(Data,Labels,2,'ploton');

%try it online
thresh = Fits{end}.ThresholdSep
[Class,Probability] = RelativeRBFOnline(Data,Fits,thresh);

Correct = Class == Labels;
Correct(Class == 0) = [];

accuracy = mean(Correct)

[ClassG,AccuracyG,ProbabilityG] = GaussianClassification(Data,Labels);
AccuracyG


%% test 2D data

means{1} = [2,2];
sigmas{1} = [1,0;0,1];

means{2} = means{1} + [dist/sqrt(2), dist/sqrt(2) ];
sigmas{2} = [1,0;0,1];

[Data,Labels] = CreateData(means,sigmas,nn);

figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
title('initial data')

[Fits] = RelativeRBFTrain(Data,Labels,2,'ploton');

%try it online
thresh = Fits{end}.ThresholdSep
[Class,Probability] = RelativeRBFOnline(Data,Fits,thresh);

Correct = Class == Labels;
Correct(Class == 0) = [];

accuracy = mean(Correct)

%try just gaussian
[ClassG,AccuracyG,ProbabilityG] = GaussianClassification(Data,Labels);
AccuracyG

%Try LDA
MdlLinear = fitcdiscr(Data,Labels);
classLDA = predict(MdlLinear,Data);
CorrLDA = classLDA == Labels;
AccuracyLDA = mean(CorrLDA)



%% test 3D data

means{1} = [2,2,2];
sigmas{1} = [1,0,0;0,1,0;0,0,1];

means{2} = means{1} + [dist/sqrt(3), dist/sqrt(3), dist/sqrt(3)];
sigmas{2} = [1,0,0;0,1,0;0,0,1];


[Data,Labels] = CreateData(means,sigmas,nn);


figure
gscatter3(Data(:,1),Data(:,2),Data(:,3),Labels,'br');
title('initial data')

[Fits] = RelativeRBFTrain(Data,Labels,2,'ploton');
 
 
 %try it online
thresh = Fits{end}.ThresholdSep
[Class,Probability] = RelativeRBFOnline(Data,Fits,thresh);

Correct = Class == Labels;
Correct(Class == 0) = [];

accuracy = mean(Correct)

%Try just gaussian classify
[ClassG,AccuracyG,ProbabilityG] = GaussianClassification(Data,Labels);
AccuracyG

%Try LDA
MdlLinear = fitcdiscr(Data,Labels);
classLDA = predict(MdlLinear,Data);
CorrLDA = classLDA == Labels;
AccuracyLDA = mean(CorrLDA)

%% test 1-D multi modal


nn = [1000,1000];

means{1} = [2];
sigmas{1} = [1];
means{2} = [6];
sigmas{2} = [1];

[Data1,Labels1] = CreateData(means,sigmas,nn);


means{1} = [12];
sigmas{1} = [1];
means{2} = [-6];
sigmas{2} = [1];

[Data2,Labels2] = CreateData(means,sigmas,nn);

Data = [Data1;Data2];
Labels = [Labels1;Labels2];

plotheight = [ones(nn(1),1)*1; ones(nn(1),1)*1.2;ones(nn(1),1)*1;ones(nn(1),1)*1.2];
figure
gscatter(Data,plotheight,Labels,'br');
title('initial data')
ylim([0,10])

[Fits] = RelativeRBFTrain(Data,Labels,2,'ploton');

%try it online
thresh = Fits{end}.ThresholdSep;
[Class,Probability] = RelativeRBFOnline(Data,Fits,thresh);

Correct = Class == Labels;
Correct(Class == 0) = [];

accuracy = mean(Correct)




%% Phase portrait data

paramNoise = 0.1; %0.2, 0.7
statenoise = 0.1;
runs = 100;
nn = 100;
x = linspace(0,2.5,nn)';

paramsinit(1,:) = [-2,1.7];
paramsinit(2,:) = [-2.4,2.0];

Data = [];
Labels = [];
for rr = 1:runs
    xon = x + randn(nn,1)*statenoise;
    
    %parameter noise
    for jj = 1:size(paramsinit,1)
        for ii = 1:length(paramsinit(jj,:))
            params(jj,ii) = paramsinit(jj,ii) + (rand(1)-.5)*paramNoise*max(paramsinit(jj,:));
        end
    end

    %Nonlinear equation
    xdot1 = -params(1,1) + params(1,1).*exp(-params(1,2).*xon);
    xdot2 = -params(2,1) + params(2,1).*exp(-params(2,2).*xon);
    
    %state noise
    xdot1 = xdot1 + randn(nn,1)*statenoise;
    xdot2 = xdot2 + randn(nn,1)*statenoise;

    Data = [Data; xon,xdot1; xon,xdot2];
    Labels = [Labels; ones(nn,1)*1; ones(nn,1)*2]; 
end


figure
gscatter(Data(:,1),Data(:,2),Labels,'br');
title('initial data')


%Train relative RBF
[Fits] = RelativeRBFTrain(Data,Labels,2,'ploton');


%try it online
thresh = Fits{end}.ThresholdSep
[Class,Probability] = RelativeRBFOnline(Data,Fits,thresh);

Correct = Class == Labels;
Correct(Class == 0) = [];
accuracy = mean(Correct)


%or try with just least squares
modelfunc = @(X) [X(:,1).^0,X(:,1),X(:,1).^2,X(:,1).^3];
outputcolumn = 2;
[LS_Param] = LeastSquaresTrain(Data,Labels,modelfunc,outputcolumn);
[ClassLS,AccuracyLS,minErrorLS] = LeastSquaresClassification(Data,Labels,LS_Param,modelfunc,outputcolumn);
AccuracyLS


%plot ls
LeastSquaresPlot(Data,Labels,LS_Param,modelfunc,[1,2]);

%% Test linear gaussian from 

modelfunc = @(X) [X(:,1).^0,X(:,1),X(:,1).^2,X(:,1).^3];
outputcolumn = 2;

[LGModel] = RBFLinearGaussian(Fits,modelfunc,outputcolumn);

%% Test weighted least squares using seperability measure

modelfunc = @(X) [X(:,1).^0,X(:,1),X(:,1).^2];
outputcolumn = 2;

[RLSmodel] = RelativeRBFWeightedLS(Data,Labels,modelfunc,outputcolumn);

[ClassRLS,AccuracyRLS,minErrorRLS] = LeastSquaresClassification(Data,Labels,RLSmodel,modelfunc,outputcolumn);
AccuracyRLS

%plot ls
LeastSquaresPlot(Data,Labels,RLSmodel,modelfunc,[1,2]);



