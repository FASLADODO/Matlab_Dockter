%With nonlinearfunctions

clear all


%%
nn = 1000;
noiz = 0.2;

rangen = [0 6*pi];

paramsnl1 = [10 1 4 0.5]; %m, b
paramsnl2 = [10 1 5 0.5]; %m, b


xm = linspace(rangen(1),rangen(2),nn)';

xm1 = xm + randn(nn,1)*noiz;
xm2 = xm + randn(nn,1)*noiz;

yx1 = ones(nn,1)*paramsnl1(1) + xm1.*paramsnl1(2) + paramsnl1(3).*sin(paramsnl1(4).*xm1) + randn(nn,1)*noiz;
yx2 = ones(nn,1)*paramsnl2(1) + xm2.*paramsnl2(2) + paramsnl2(3).*cos(paramsnl2(4).*xm2) + randn(nn,1)*noiz;


X1 = [xm2, yx1];
X2 = [xm1, yx2];

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


std(X1)
std(X2)

(range(Data)./length(Data)) .* std(Data)

mean(sqrt(std(Data)))

%% plot the resulting seperations

[h] = windowDiscriminantPlot(Data, Y, 0.5);


%% Train the models for each
pthresh = 0.5;
kmeanlimit = 7;

if(~exist('kmg'))
    kmg = [];
end
[Models, kmg] = windowDiscriminantTrain(Data, Y, pthresh, kmg, kmeanlimit, 'ploton');


%% Now get some new online data

rng('shuffle')

xm1_on = xm + randn(nn,1)*noiz;
xm2_on = xm + randn(nn,1)*noiz;

yx1_on = ones(nn,1)*paramsnl1(1) + xm1_on.*paramsnl1(2) + paramsnl1(3).*sin(paramsnl1(4).*xm1_on) + randn(nn,1)*noiz;
yx2_on = ones(nn,1)*paramsnl2(1) + xm2_on.*paramsnl2(2) + paramsnl2(3).*cos(paramsnl2(4).*xm2_on) + randn(nn,1)*noiz;


X1_Online = [xm2_on, yx1_on];
X2_Online = [xm1_on, yx2_on];

Data_Online = [X1_Online;X2_Online];


%% Log Likelihood ratios

disp('classification accuracy with Rods Discriminant:')
LLR_X1 = windowDiscriminantOnline(X1_Online, Models, classes, pthresh);
LLR_X1(end);
class_okay1 = LLR_X1 >= 0;
Accuracy1 = sum(class_okay1)/length(class_okay1)

LLR_X2 = windowDiscriminantOnline(X2_Online, Models, classes, pthresh);
LLR_X2(end);
class_okay2 = LLR_X2 <= 0;
Accuracy2 = sum(class_okay2)/length(class_okay2)

disp('classification accuracy with gaussians:')
AccLC = LeastSqauresClassify( [Data_Online(:,1),ones(length(Data_Online),1)], Data_Online(:,2), Y);

figure
plot(LLR_X1,'r');
figure
plot(LLR_X2,'b');











    